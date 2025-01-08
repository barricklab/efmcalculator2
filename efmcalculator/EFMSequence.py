from Bio.SeqRecord import SeqRecord


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

        self.formatted_name = None
        self.is_circular = is_circular

        self._ssrs = None
        self._srss = None
        self._rmds = None

        self._unique_annotations = {}

        self._selected_predictions = []
        self._selected_annotations = []
        self._filtered_ssrs = None
        self._filtered_srss = None
        self._filtered_rmds = None

        self._originhash = originhash

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
        if self._unique_annotations == {}:
            features = sequence_to_features_df(self.seq, self.circular)
        return self._unique_annotations

    def call_predictions(self, strategy):
        seq = str(self.seq.strip("\n").upper().replace("U", "T"))
        ssr_df, srs_df, rmd_df = predict(seq, strategy, self.is_circular)
        self._ssrs, self._srss, self._rmds = post_process(ssr_df,
                                                             srs_df,
                                                             rmd_df,
                                                             self,
                                                             self.is_circular)

    def same_origin(self, other):
        if not self._originhash or not other._filehash:
            return False
        elif self._originhash != other._filehash:
            return False
        else:
            return True
