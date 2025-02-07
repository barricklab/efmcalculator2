from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SimpleLocation, CompoundLocation, SeqFeature

import polars as pl
from ..webapp.vis_utils import eval_top
import streamlit as st

from ..pipeline.features import seqfeature_hash, sequence_to_features_df

from ..pipeline.primary_pipeline import predict
from ..pipeline.post_process import post_process
from ..pipeline.features import sequence_to_features_df

class EFMSequence(SeqRecord):
    """SeqRecord child class with equality handling and prediction methods"""
    def __init__(self, seqrecord: SeqRecord, is_circular = False, originhash = None):
        super().__init__(seq=seqrecord.seq,
            id=seqrecord.id,
            name=seqrecord.name,
            description=seqrecord.description,
            dbxrefs = seqrecord.dbxrefs,
            features = seqrecord.features,
            annotations = seqrecord.annotations,
            letter_annotations = seqrecord.letter_annotations)
        self.seq = seqrecord.seq
        self.formatted_name = None
        self.is_circular = is_circular

        self._predicted = False
        self._top = None
        self._ssrs = None
        self._srss = None
        self._rmds = None

        self._unique_annotations = {}
        self._originhash = originhash

        self._filtered_ssrs = None
        self._filtered_srss = None
        self._filtered_rmds = None
        self._filtered_top = None
        self._last_filters = []

        self._plotted_predictions = []


    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, seq):
        self._ssrs = None
        self._srss = None
        self._rmds = None
        self._predicted = False
        self._seq = seq

    @property
    def predicted(self):
        return self._predicted

    @property
    def top(self):
        if self._top is None:
            raise ValueError("No predictions calculated. Run call_predictions() first.")
        return self._top

    @property
    def ssrs(self):
        if self._ssrs is None:
            raise ValueError("No SSRs calculated. Run call_predictions() first.")
        return self._ssrs

    @property
    def srss(self):
        if self._srss is None:
            raise ValueError("No SRSs calculated. Run call_predictions() first.")
        return self._srss

    @property
    def rmds(self):
        if self._rmds is None:
            raise ValueError("No RMDs calculated. Run call_predictions() first.")
        return self._rmds

    @property
    def unique_annotations(self):
        if not self.annotations:
            return None
        if isinstance(self._unique_annotations, dict) and self._unique_annotations == {}:
            self._unique_annotations = sequence_to_features_df(self, self.is_circular)
            self._unique_annotations = self._unique_annotations.with_columns(
                pl.concat_str([pl.col("annotations"), pl.lit(" ("), pl.col("left_bound"), pl.lit("-"), pl.col("right_bound"), pl.lit(")")]).alias("annotationobjexpanded_names")
            )
        unique_expaned_names = self._unique_annotations.select(pl.col("annotationobjexpanded_names")).unique().rows()
        unique_expaned_names = [x[0] for x in unique_expaned_names]

        return sorted(unique_expaned_names)

    def call_predictions(self, strategy):
        seq = str(self.seq.strip("\n").upper().replace("U", "T"))
        ssr_df, srs_df, rmd_df = predict(seq, strategy, self.is_circular)
        ssr_df, srs_df, rmd_df = post_process(ssr_df,
                                                             srs_df,
                                                             rmd_df,
                                                             self,
                                                             self.is_circular)
        self._ssrs = ssr_df.with_columns(
            pl.concat_str(pl.col(['repeat', 'repeat_len', "start", "count", "mutation_rate"])).hash().cast(pl.String).alias('predid')
        )
        self._srss = srs_df.with_columns(
            pl.concat_str(pl.col(['repeat', 'repeat_len', "first_repeat", "second_repeat", "distance", "mutation_rate"])).hash().cast(pl.String).alias('predid')
        )
        self._rmds = rmd_df.with_columns(
            pl.concat_str(pl.col(['repeat', 'repeat_len', "first_repeat", "second_repeat", "distance", "mutation_rate"])).hash().cast(pl.String).alias('predid')
        )
        self._top = eval_top(self._ssrs, self._srss, self._rmds)
        self._filtered_top = self._top
        self._plotted_predictions = [x[0] for x in self._top.select(pl.col("predid")).unique().rows()]
        self.set_filters([])
        self._predicted = True

    def set_filters(self, annotations):
        filters_changed = False
        if self._last_filters != annotations:
            filters_changed = True
            self._last_filters = annotations

        if annotations:
            annotation_objects = self._unique_annotations.filter(pl.col("annotationobjexpanded_names").is_in(annotations))
            annotation_objects = annotation_objects.select(pl.col("annotationobjects")).unique().rows()
            annotation_objects = [x[0] for x in annotation_objects]
            self._filtered_ssrs = self._ssrs.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
            self._filtered_srss = self._srss.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
            self._filtered_rmds = self._rmds.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
            self._filtered_top = eval_top(self._filtered_ssrs, self._filtered_srss, self._filtered_rmds)
        else:
            self._filtered_ssrs = self._ssrs
            self._filtered_srss = self._srss
            self._filtered_rmds = self._rmds
            self._filtered_top = self._top

        self._filtered_ssrs = self._filtered_ssrs.with_columns(
            pl.when(pl.col("predid").is_in(self._plotted_predictions))
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("show")
        ).sort(by="mutation_rate", descending=True)
        self._filtered_srss = self._filtered_srss.with_columns(
            pl.when(pl.col("predid").is_in(self._plotted_predictions))
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("show")
        ).sort(by="mutation_rate", descending=True)
        self._filtered_rmds = self._filtered_rmds.with_columns(
            pl.when(pl.col("predid").is_in(self._plotted_predictions))
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("show")
        ).sort(by="mutation_rate", descending=True)
        self._filtered_top = self._filtered_top.with_columns(
            pl.when(pl.col("predid").is_in(self._plotted_predictions))
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("show")
        ).sort(by="mutation_rate", descending=True)

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

    def upate_top_session(self):
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
    def update_ssr_session(self):
        changes = st.session_state["ssrchanges"]['edited_rows']
        for change in changes:
            try:
                new_state = changes[change]['show']
                changed_id = self._filtered_ssrs[change]['predid'][0]
                if new_state:
                    self._plotted_predictions.append(changed_id)
                else:
                    self._plotted_predictions.remove(changed_id)
            except ValueError:
                pass
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

    def same_origin(self, other):
        if not self._originhash or not other._filehash:
            return False
        elif self._originhash != other._filehash:
            return False
        else:
            return True

def sequence_to_features_df(sequence, circular=True):
    """Takes genbank annotations and turns them into a polars dataframe"""
    features = sequence.features
    seqlen = len(sequence)

    def get_feature_bounds(feature_location):
        if isinstance(feature_location, SimpleLocation):
            return feature_location.start, feature_location.end
        elif isinstance(feature_location, CompoundLocation):
            if feature_location.parts[-1].end < feature_location.parts[0].start:
                end = len(sequence) + feature_location.parts[-1].end
            else:
                end = feature_location.parts[-1].end
            return feature_location.parts[0].start, end

    if not circular:
        # Look for compound wraparounds and break them into simple features
        newfeatures = []
        deletedfeatures = []
        for feature in features:
            if not isinstance(feature.location, CompoundLocation):
                continue

            # Check whether the compound feature is actually a wraparound
            wraparound_part_index = None
            last_part_start = None
            rightmost_part = None
            for i, part in enumerate(feature.location.parts):
                if rightmost_part != None and part.start < last_part_start:
                    wraparound_part_index = i
                if rightmost_part == None:
                    rightmost_part = i
                    last_part_start = part.start
            if wraparound_part_index is None:
                continue

            # If it is a wraparound, break it into two features
            leftsplit_locations = feature.location.parts[:wraparound_part_index]
            if len(leftsplit_locations) > 1:
                leftsplit_locations = CompoundLocation(leftsplit_locations)
            else:
                leftsplit_locations = leftsplit_locations[0]
            newleftfeature = SeqFeature(leftsplit_locations, feature.type, qualifiers=feature.qualifiers)


            rightsplit_locations = feature.location.parts[wraparound_part_index:]
            if len(rightsplit_locations) > 1:
                rightsplit_locations = CompoundLocation(rightsplit_locations)
            else:
                rightsplit_locations = rightsplit_locations[0]
            newrightfeature = SeqFeature(rightsplit_locations, feature.type, qualifiers=feature.qualifiers)

            newfeatures.append(newleftfeature)
            newfeatures.append(newrightfeature)
            deletedfeatures.append(feature)

        for feature in deletedfeatures:
            features.remove(feature)
        features.extend(newfeatures)

    df = pl.DataFrame([(feature.type,
                        get_feature_bounds(feature.location),
                        feature.qualifiers.get("label", "")[0],
                        seqfeature_hash(feature)) for feature in features],
        schema=['type', 'loc', 'annotations', 'annotationobjects'],
        orient="row")

    # expand out loc
    df = df.with_columns(pl.col("loc").list.to_struct(fields=['left_bound', 'right_bound'])).unnest("loc")
    df = df.with_columns(
        pl.col("left_bound").cast(pl.Int32),
        pl.col("right_bound").cast(pl.Int32)
    )

    return df
