from Bio.SeqRecord import SeqRecord

import polars as pl
from .short_seq_finder import predict
from .post_process import post_process
from .features import sequence_to_features_df

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
        self._ssrs, self._srss, self._rmds = post_process(ssr_df,
                                                             srs_df,
                                                             rmd_df,
                                                             self,
                                                             self.is_circular)
        self._predicted = True


    def set_filters(self, annotations):
        if annotations:
            annotation_objects = self._unique_annotations.filter(pl.col("annotationobjexpanded_names").is_in(annotations))
            annotation_objects = annotation_objects.select(pl.col("annotationobjects")).unique().rows()
            annotation_objects = [x[0] for x in annotation_objects]
            self._filtered_ssrs = self._ssrs.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
            self._filtered_srss = self._srss.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
            self._filtered_rmds = self._rmds.filter(pl.col("annotationobjects").list.set_intersection(annotation_objects).list.len() != 0)
        else:
            self._filtered_ssrs = self._ssrs
            self._filtered_srss = self._srss
            self._filtered_rmds = self._rmds


    def update_ssr_session(self, changes, from_top=False):
        pass
    def update_srs_session(self, changes, from_top=False):
        pass
    def update_rmd_session(self, changes, from_top=False):
        pass

    def same_origin(self, other):
        if not self._originhash or not other._filehash:
            return False
        elif self._originhash != other._filehash:
            return False
        else:
            return True
