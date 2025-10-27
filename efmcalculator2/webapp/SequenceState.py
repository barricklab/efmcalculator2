"""Manages sequence-specific state information for the webapp"""
import polars as pl
from .vis_utils import eval_top
import streamlit as st
import logging
for name, l in logging.root.manager.loggerDict.items():
    if "streamlit" in name:
        l.disabled = True

from ..ingest.EFMSequence import sequence_to_features_df
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode, JsCode

class SequenceState():
    """Manages sequence-specific state information for the webapp"""
    def __init__(self, efmsequence):
        self.efmsequence = efmsequence


        self._ssr_webapp_state = None

        self._webapp_ssrs = None
        self._webapp_srss = None
        self._webapp_rmds = None
        self._webapp_top = None

        self._unique_annotations = {}

        self.ready = False
        self.refreshed = False
        self.counter = 0

        self.topkey = f"toptable{self.counter}"
        self.ssrkey = f"ssrtable{self.counter}"
        self.srskey = f"srstable{self.counter}"
        self.rmdkey = f"rmdtable{self.counter}"

    # Data from upstream EFMSequence

    @property
    def predicted(self):
        return self.efmsequence.predicted

    @property
    def unique_annotations(self):
        return self.efmsequence.unique_annotations

    @property
    def seq(self):
        return self.efmsequence.seq

    @property
    def features(self):
        return self.efmsequence.features

    def post_predict_processing(self):
        if not self.predicted:
            raise ValueError("Predictions have not yet been completed")
        elif self.ready:
            return
        self._webapp_ssrs = self.efmsequence._ssrs.sort(by="mutation_rate", descending=True)
        self._webapp_srss = self.efmsequence._srss.sort(by="mutation_rate", descending=True)
        self._webapp_rmds = self.efmsequence._rmds.sort(by="mutation_rate", descending=True)
        self._webapp_top = self.efmsequence._top.sort(by="mutation_rate", descending=True)
        self._shown_predictions = [x[0] for x in self.efmsequence._top.select(pl.col("predid")).unique().rows()]
        self._filtered_ssrs = None
        self._filtered_srss = None
        self._filtered_rmds = None
        self._filtered_top = None
        self._last_filters = []
        self.set_filters(None)

        self.last_top_selections = None
        self.last_ssr_selections = None
        self.last_srs_selections = None
        self.last_rmd_selections = None

    def reset_selected_predictions(self):
        self._shown_predictions = [x[0] for x in self._filtered_top.select(pl.col("predid")).unique().rows()]

    @property
    def unique_annotations(self):
        if not self.efmsequence.annotations:
            return None
        if isinstance(self._unique_annotations, dict) and self._unique_annotations == {}:
            self._unique_annotations = sequence_to_features_df(self.efmsequence, self.efmsequence.is_circular)
            self._unique_annotations = self._unique_annotations.with_columns(
                pl.concat_str([pl.col("annotations"), pl.lit(" ("), pl.col("left_bound"), pl.lit("-"), pl.col("right_bound"), pl.lit(")")]).alias("annotationobjexpanded_names")
            )
        unique_expaned_names = self._unique_annotations.select(pl.col("annotationobjexpanded_names")).unique().rows()
        unique_expaned_names = [x[0] for x in unique_expaned_names]

        return sorted(unique_expaned_names)

    # Data pertainint to webapp
    def set_filters(self, annotations):
        if annotations:
            annotation_objects = self._unique_annotations.filter(pl.col("annotationobjexpanded_names").is_in(annotations))
            annotation_objects = annotation_objects.select(pl.col("annotationobjects")).unique().rows()
            annotation_objects = [x[0] for x in annotation_objects]
            self._filtered_ssrs = self._webapp_ssrs.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
            self._filtered_srss = self._webapp_srss.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
            self._filtered_rmds = self._webapp_rmds.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
            self._filtered_top = eval_top(self._filtered_ssrs, self._filtered_srss, self._filtered_rmds).with_columns(pl.lit(False).alias("show"))
        else:
            self._filtered_ssrs = self._webapp_ssrs
            self._filtered_srss = self._webapp_srss
            self._filtered_rmds = self._webapp_rmds
            self._filtered_top = self._webapp_top

    @property
    def top_webapp_table(self):
        """Streamlit aggrid table representing ssr data"""

        pandas_conversion = self._filtered_top.to_pandas()
        return AgGrid(pandas_conversion,
                            gridOptions=self._top_webapp_state,
                            height=500,
                            fit_columns_on_grid_load=True,
                            allow_unsafe_jscode = True,
                            update_mode = GridUpdateMode.SELECTION_CHANGED,
                            key = self.topkey,
                            callback = self.update_top_session)


    def rebuild_top_table(self):
        cell_hover_handler = JsCode("""""")
        js_hover_handler = """"""

        preselected_indices = self.get_shown_predictions(self._filtered_top)
        preselected_indices = [str(index) for index in preselected_indices]
        if self.last_top_selections is None:
            self.last_top_selections = preselected_indices

        pandas_conversion = self._filtered_top.to_pandas()
        builder = GridOptionsBuilder.from_dataframe(pandas_conversion)
        builder.configure_selection(selection_mode='multiple', use_checkbox= True, pre_selected_rows= preselected_indices)

        builder.configure_grid_options(onCellMouseOver=cell_hover_handler)
        builder.configure_column("repeat", header_name="Sequence", tooltipField="repeat")
        builder.configure_column("mutation_rate", header_name="Mutation Rate",
                    type=["numericColumn"], valueFormatter="x.toExponential(2)")
        builder.configure_column("annotations", header_name="Annotations", tooltipField="annotations")
        builder.configure_column("predid", hide = True)
        builder.configure_column("annotationobjects", hide = True)
        grid_options = builder.build()
        self._top_webapp_state = grid_options

    def update_top_session(self, callbackobj):
        self.refresh_last_shown()
        state = callbackobj.grid_response["gridState"]
        new_selection = state.get('rowSelection', [])

        if self.last_top_selections is not None and self.last_top_selections == new_selection:
            return

        if self.last_top_selections is None:
            self.last_top_selections = []
        if self.last_top_selections == new_selection:
            return

        self.counter += 1
        self.ssrkey = f"ssrtable{self.counter}"
        self.srskey = f"srstable{self.counter}"
        self.rmdkey = f"rnmdable{self.counter}"

        to_drop = [int(x) for x in self.last_top_selections if x not in new_selection]
        to_add = [int(x) for x in new_selection if x not in self.last_top_selections]


        self._update_general_table(self._filtered_top, to_drop, to_add)

        self.last_top_selections = new_selection


    # SSR Table Functions
    @property
    def ssr_webapp_table(self):
        """Streamlit aggrid table representing ssr data"""

        pandas_conversion = self._filtered_ssrs.to_pandas()
        return AgGrid(pandas_conversion,
                            gridOptions=self._ssr_webapp_state,
                            height=500,
                            fit_columns_on_grid_load=True,
                            allow_unsafe_jscode = True,
                            update_mode = GridUpdateMode.SELECTION_CHANGED,
                            key=self.ssrkey,
                            callback=self.update_ssr_session,
                            reload_data=False)

    def rebuild_ssr_table(self):
        cell_hover_handler = JsCode("""
                    function(params) {
                        // debug
                        const clickedColumn = params.column.colId;
                        const clickedRowIndex = params.rowIndex;
                        const clickedValue = params.node.data[clickedColumn];

                        const predidValue = params.node.data["predid"];

                        // Display information about the click
                        const message = `You hovered on row ${clickedRowIndex}, column ${clickedColumn}, value is ${predidValue}`;
                        console.log(message);

                        window.parent.postMessage(
                            {type: 'cellMouseOver', predid: predidValue}, '*'
                            );

                    }
                    """)
        js_hover_handler = """
                        window.addEventListener("message", (event) => {
                            if (event.data.type === "cellMouseOver") {
                                resolve(event.data.predid);
                            }
                        });
                    """

        preselected_indices = self.get_shown_predictions(self._filtered_ssrs)
        preselected_indices = [str(index) for index in preselected_indices]
        if self.last_ssr_selections is None:
            self.last_ssr_selections = preselected_indices

        pandas_conversion = self._filtered_ssrs.to_pandas()
        builder = GridOptionsBuilder.from_dataframe(pandas_conversion)
        builder.configure_selection(selection_mode='multiple', use_checkbox= True, pre_selected_rows= preselected_indices)

        builder.configure_grid_options(onCellMouseOver=cell_hover_handler)
        builder.configure_column("repeat", header_name="Sequence", tooltipField="repeat")
        builder.configure_column("repeat_len", header_name="Repeat Length", type=["numericColumn"])
        builder.configure_column("start", header_name="Start", type=["numericColumn"])
        builder.configure_column("count", header_name="Count", type=["numericColumn"])
        builder.configure_column("mutation_rate", header_name="Mutation Rate",
                    type=["numericColumn"], valueFormatter="x.toExponential(2)")
        builder.configure_column("annotations", header_name="Annotations", tooltipField="annotations")
        builder.configure_column("predid", hide = True)
        builder.configure_column("annotationobjects", hide = True)
        grid_options = builder.build()
        self._ssr_webapp_state = grid_options

    def update_ssr_session(self, callbackobj):
        self.refresh_last_shown()
        state = callbackobj.grid_response["gridState"]
        new_selection = state.get('rowSelection', [])

        if self.last_ssr_selections is not None and self.last_ssr_selections == new_selection:
            return
        if self.last_ssr_selections is None:
            self.last_ssr_selections = []
        if self.last_ssr_selections == new_selection:
            return

        self.counter += 1
        self.topkey = f"toptable{self.counter}"

        to_drop = [int(x) for x in self.last_ssr_selections if x not in new_selection]
        to_add = [int(x) for x in new_selection if x not in self.last_ssr_selections]


        self._update_general_table(self._filtered_ssrs, to_drop, to_add)

        self.last_ssr_selections = new_selection

    # SRS Table Functions
    @property
    def srs_webapp_table(self):
        """Streamlit aggrid table representing srs data"""

        pandas_conversion = self._filtered_srss.to_pandas()
        return AgGrid(pandas_conversion,
                            gridOptions=self._srs_webapp_state,
                            height=500,
                            fit_columns_on_grid_load=True,
                            allow_unsafe_jscode = True,
                            update_mode = GridUpdateMode.SELECTION_CHANGED,
                            key=self.srskey,
                            callback=self.update_srs_session,
                            reload_data=False)

    def rebuild_srs_table(self):
        cell_hover_handler = JsCode("""
                    function(params) {
                        // debug
                        const clickedColumn = params.column.colId;
                        const clickedRowIndex = params.rowIndex;
                        const clickedValue = params.node.data[clickedColumn];

                        const predidValue = params.node.data["predid"];

                        // Display information about the click
                        const message = `You hovered on row ${clickedRowIndex}, column ${clickedColumn}, value is ${predidValue}`;
                        console.log(message);

                        window.parent.postMessage(
                            {type: 'cellMouseOver', predid: predidValue}, '*'
                            );

                    }
                    """)
        js_hover_handler = """

                    """

        preselected_indices = self.get_shown_predictions(self._filtered_srss)
        preselected_indices = [str(index) for index in preselected_indices]
        if self.last_srs_selections is None:
            self.last_srs_selections = preselected_indices

        pandas_conversion = self._filtered_srss.to_pandas()
        builder = GridOptionsBuilder.from_dataframe(pandas_conversion)
        builder.configure_selection(selection_mode='multiple', use_checkbox= True, pre_selected_rows= preselected_indices)
        builder.configure_grid_options(onCellMouseOver=cell_hover_handler)
        builder.configure_column("repeat", header_name="Sequence", tooltipField="repeat")
        builder.configure_column("repeat_len", header_name="Repeat Length", type=["numericColumn"])
        builder.configure_column("first_repeat", header_name="First Repeat", type=["numericColumn"])
        builder.configure_column("second_repeat", header_name="Second Repeat", type=["numericColumn"])
        builder.configure_column("distance", header_name="Distance", type=["numericColumn"])
        builder.configure_column("mutation_rate", header_name="Mutation Rate",
                    type=["numericColumn"], valueFormatter="x.toExponential(2)")
        builder.configure_column("annotations", header_name="Annotations", tooltipField="annotations")
        builder.configure_column("predid", hide = True)
        builder.configure_column("annotationobjects", hide = True)

        grid_options = builder.build()
        self._srs_webapp_state = grid_options

    def update_srs_session(self, callbackobj):
        self.refresh_last_shown()
        state = callbackobj.grid_response["gridState"]
        new_selection = state.get('rowSelection', [])

        if self.last_srs_selections is not None and self.last_srs_selections == new_selection:
            return
        if self.last_srs_selections is None:
            self.last_srs_selections = []
        if self.last_srs_selections == new_selection:
            return

        self.counter += 1
        self.topkey = f"toptable{self.counter}"

        to_drop = [int(x) for x in self.last_srs_selections if x not in new_selection]
        to_add = [int(x) for x in new_selection if x not in self.last_srs_selections]


        self._update_general_table(self._filtered_srss, to_drop, to_add)

        self.last_srs_selections = new_selection

    # RMD Table Functions
    @property
    def rmd_webapp_table(self):
        """Streamlit aggrid table representing rmd data"""

        pandas_conversion = self._filtered_rmds.to_pandas()
        return AgGrid(pandas_conversion,
                            gridOptions=self._rmd_webapp_state,
                            height=500,
                            fit_columns_on_grid_load=True,
                            allow_unsafe_jscode = True,
                            update_mode = GridUpdateMode.SELECTION_CHANGED,
                            key=self.rmdkey,
                            callback=self.update_rmd_session,
                            reload_data=False)

    def rebuild_rmd_table(self):
        cell_hover_handler = JsCode("""
                    function(params) {
                        // debug
                        const clickedColumn = params.column.colId;
                        const clickedRowIndex = params.rowIndex;
                        const clickedValue = params.node.data[clickedColumn];

                        const predidValue = params.node.data["predid"];

                        // Display information about the click
                        const message = `You hovered on row ${clickedRowIndex}, column ${clickedColumn}, value is ${predidValue}`;
                        console.log(message);

                        window.parent.postMessage(
                            {type: 'cellMouseOver', predid: predidValue}, '*'
                            );

                    }
                    """)
        js_hover_handler = """

                    """

        preselected_indices = self.get_shown_predictions(self._filtered_rmds)
        preselected_indices = [str(index) for index in preselected_indices]
        if self.last_rmd_selections is None:
            self.last_rmd_selections = preselected_indices

        pandas_conversion = self._filtered_rmds.to_pandas()
        builder = GridOptionsBuilder.from_dataframe(pandas_conversion)
        builder.configure_selection(selection_mode='multiple', use_checkbox= True, pre_selected_rows= preselected_indices)
        builder.configure_grid_options(onCellMouseOver=cell_hover_handler)
        builder.configure_column("repeat", header_name="Sequence", tooltipField="repeat")
        builder.configure_column("repeat_len", header_name="Repeat Length", type=["numericColumn"])
        builder.configure_column("first_repeat", header_name="First Repeat", type=["numericColumn"])
        builder.configure_column("second_repeat", header_name="Second Repeat", type=["numericColumn"])
        builder.configure_column("distance", header_name="Distance", type=["numericColumn"])
        builder.configure_column("mutation_rate", header_name="Mutation Rate",
                    type=["numericColumn"], valueFormatter="x.toExponential(2)")
        builder.configure_column("annotations", header_name="Annotations", tooltipField="annotations")
        builder.configure_column("predid", hide = True)
        builder.configure_column("annotationobjects", hide = True)

        grid_options = builder.build()
        self._rmd_webapp_state = grid_options

    def update_rmd_session(self, callbackobj):
        self.refresh_last_shown()
        state = callbackobj.grid_response["gridState"]
        new_selection = state.get('rowSelection', [])

        if self.last_rmd_selections is not None and self.last_rmd_selections == new_selection:
            return
        if self.last_rmd_selections is None:
            self.last_rmd_selections = []
        if self.last_rmd_selections == new_selection:
            return

        self.counter += 1
        self.topkey = f"toptable{self.counter}"

        to_drop = [int(x) for x in self.last_rmd_selections if x not in new_selection]
        to_add = [int(x) for x in new_selection if x not in self.last_rmd_selections]


        self._update_general_table(self._filtered_rmds, to_drop, to_add)

        self.last_rmd_selections = new_selection


    def annotation_coverage(self, annotations):
        annotation_objects = self._unique_annotations.filter(pl.col("annotationobjexpanded_names").is_in(annotations))
        annotation_objects = annotation_objects.select(["left_bound", "right_bound"])

        coverage = []
        for row in annotation_objects.iter_rows(named=True):
            for i, occupied_area in enumerate(coverage):
                if occupied_area[0] <= row['left_bound'] and occupied_area[1] >= row['right_bound']:
                    # Entirely inside
                    break
                elif occupied_area[0] <= row['left_bound'] and row['left_bound'] <= occupied_area[1] <= row['right_bound']:
                    coverage[i][0] = row['left_bound']
                    break
                elif  row['left_bound'] <= occupied_area[0] <= row['right_bound'] and occupied_area[1] >= row['right_bound']:
                    coverage[i][1] = row['right_bound']
                    break
            else:
                # entirely outside
                coverage.append((row['left_bound'], row['right_bound']))
        base_coverage = 0
        for region in coverage:
            base_coverage += region[1] - region[0] + 1
        return base_coverage


    def get_shown_predictions(self, df):
        selected_rows = df.with_row_index().with_columns(
            pl.when(pl.col("predid").is_in(self._shown_predictions))
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("show")
        ).filter(pl.col("show") == True).select("index").unique().rows()
        selected_rows = [x for xs in selected_rows for x in xs] # Unnest
        return selected_rows

    def _update_general_table(self, dataframe, to_drop, to_add):
        """Internal function. Updates shown predictions based on those selected in table"""
        to_drop = dataframe.with_row_index(
                    ).filter(pl.col("index").is_in(to_drop)
                    ).select('predid').rows()
        to_drop = [x for xs in to_drop for x in xs] # Unnest

        to_add = dataframe.with_row_index(
                    ).filter(pl.col("index").is_in(to_add)
                    ).select('predid').rows()
        to_add = [x for xs in to_add for x in xs] # Unnest

        for item in to_drop:
            try:
                self._shown_predictions.remove(item)
            except ValueError:
                pass # Eronious removal
        for item in to_add:
            if item not in self._shown_predictions:
                self._shown_predictions.append(item)

    def refresh_last_shown(self):
        if not self.refreshed:
            #top
            preselected_indices = self.get_shown_predictions(self._filtered_top)
            preselected_indices = [str(index) for index in preselected_indices]
            self.last_top_selections = preselected_indices
            #ssrs
            preselected_indices = self.get_shown_predictions(self._filtered_ssrs)
            preselected_indices = [str(index) for index in preselected_indices]
            self.last_ssr_selections = preselected_indices
            #srss
            preselected_indices = self.get_shown_predictions(self._filtered_srss)
            preselected_indices = [str(index) for index in preselected_indices]
            self.last_srs_selections = preselected_indices
            #rmd
            preselected_indices = self.get_shown_predictions(self._filtered_rmds)
            preselected_indices = [str(index) for index in preselected_indices]
            self.last_rmd_selections = preselected_indices
            self.refreshed = True
