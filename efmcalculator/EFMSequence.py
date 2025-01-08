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
        self.seq = seqrecord.seq
        self.formatted_name = None
        self.is_circular = is_circular

        self._predicted = False
        self._ssrs = None
        self._srss = None
        self._rmds = None

        self._unique_annotations = {}

        self.session_dataframes = []

        self._originhash = originhash

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
        self._predicted = True

    def update_ssr_session(self, ssr_df):
        self.session_dataframes[0] = ssr_df
    def update_srs_session(self, srs_df):
        self.session_dataframes[1] = srs_df
    def update_rmd_session(self, rmd_df):
        self.session_dataframes[2] = rmd_df

    def same_origin(self, other):
        if not self._originhash or not other._filehash:
            return False
        elif self._originhash != other._filehash:
            return False
        else:
            return True
