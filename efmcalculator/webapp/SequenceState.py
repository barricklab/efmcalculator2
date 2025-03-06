"""Manages sequence-specific state information for the webapp"""
import polars as pl
from .vis_utils import eval_top
import streamlit as st
from rich import print

from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode, JsCode


class SequenceState():
    """Manages sequence-specific state information for the webapp"""
    def __init__(self, efmsequence):
        self.efmsequence = efmsequence

        self._shown_predictions = [x[0] for x in self.efmsequence._top.select(pl.col("predid")).unique().rows()]
        self._filtered_ssrs = None
        self._filtered_srss = None
        self._filtered_rmds = None
        self._filtered_top = None
        self._last_filters = []
        self.set_filters(None)

        self._ssr_webapp_state = None

        self._webapp_ssrs = self.efmsequence._ssrs
        self._webapp_srss = self.efmsequence._srss
        self._webapp_rmds = self.efmsequence._rmds
        self._webapp_top = self.efmsequence._top

    # Data from upstream EFMSequence

    @property
    def predicted(self):
        return self.efmsequence.predicted

    # Data pertainint to webapp
    def set_filters(self, annotations):
        if annotations:
            annotation_objects = self.efmsequence._unique_annotations.filter(pl.col("annotationobjexpanded_names").is_in(annotations))
            annotation_objects = annotation_objects.select(pl.col("annotationobjects")).unique().rows()
            annotation_objects = [x[0] for x in annotation_objects]
            self._filtered_ssrs = self._webapp_ssrs.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
            self._filtered_srss = self._webapp_srss.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
            self._filtered_rmds = self._webapp_rmds.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
            self._filtered_top = eval_top(self._filtered_ssrs, self._filtered_srss, self._filtered_rmds).with_columns(pl.lit(False).alias("show")).sort(by="mutation_rate", descending=True)
        else:
            self._filtered_ssrs = self._webapp_ssrs
            self._filtered_srss = self._webapp_srss
            self._filtered_rmds = self._webapp_rmds
            self._filtered_top = self._webapp_top

    @property
    def ssr_webapp_table(self):
        """Streamlit aggrid table representing ssr data"""

        cell_hover_handler = JsCode("""""")
        js_hover_handler = """"""
        pandas_conversion = self._filtered_ssrs.to_pandas()
        builder = GridOptionsBuilder.from_dataframe(pandas_conversion)

        preselected_indices = pandas_conversion[pandas_conversion["show"] == True].index.tolist()
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

        return AgGrid(pandas_conversion,
                            gridOptions=grid_options,
                            height=500,
                            fit_columns_on_grid_load=True,
                            allow_unsafe_jscode = True,
                            update_mode = GridUpdateMode.SELECTION_CHANGED)

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

    def update_top_session(self):
        changes = st.session_state["topchanges"]['edited_rows']
        for change in changes:
            try:
                new_state = changes[change]['show']
                changed_id = self._filtered_top[change]['predid'][0]
                if new_state:
                    self._plotted_predictions.append(changed_id)
                else:
                    self._plotted_predictions.remove(changed_id)
            except ValueError:
                pass

    def update_ssr_session(self, options, state): #https://discuss.streamlit.io/t/persist-aggrid-state-across-page-changes/37882/3

        # Extract current state
        keys = (('aggregation', 'aggregationModel'), ('columnSizing', 'columnSizingModel'), ('sort', 'sortModel'))
        groups = state.get('rowGroup',{}).get('groupColIds',[])
        order = state.get('columnOrder',{}).get('orderedColIds',[])
        hidden = state.get('columnVisibility', {}).get('hiddenColIds', [])
        selected = state.get('rowSelection')
        fields = {}
        for key in keys:
            if a:= state.get(key[0]):
                a = a[key[1]]
                for c in a:
                    col = c.pop('colId')
                    try:
                        fields[col].update(c)
                    except KeyError:
                        fields[col] = c

        # Save state
        for c in options:
            c['sort'] = ''
            if c['field'] in fields:
                c.update(fields[c['field']])
                if order:
                    c['order'] = order.index(c['field'])
            if c['field'] in groups:
                c['rowGroup'] = True
            else:
                c['rowGroup'] = False
            if c['field'] in hidden:
                c['hide'] = True
            else:
                c['hide'] = False
        print(options)
        self._ssr_webapp_state = options

    def update_srs_session(self):
        changes = st.session_state["srschanges"]['edited_rows']
        for change in changes:
            try:
                new_state = changes[change]['show']
                changed_id = self._filtered_srss[change]['predid'][0]
                if new_state:
                    self._plotted_predictions.append(changed_id)
                else:
                    self._plotted_predictions.remove(changed_id)
            except ValueError:
                pass
    def update_rmd_session(self):
        changes = st.session_state["rmdchanges"]['edited_rows']
        for change in changes:
            try:
                new_state = changes[change]['show']
                changed_id = self._filtered_rmds[change]['predid'][0]
                if new_state:
                    self._plotted_predictions.append(changed_id)
                else:
                    self._plotted_predictions.remove(changed_id)
            except ValueError:
                pass

    def get_shown_predictions(self, df):
        selected_rows = df.with_row_index().with_columns(
            pl.when(pl.col("predid").is_in(self._shown_predictions))
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("show")
        ).filter(pl.col("show") == True).select("index").unique()

        return selected_rows
