from Bio.SeqRecord import SeqRecord


from .short_seq_finder import predict
from .post_process import post_process

class EFMSequence(SeqRecord):
    """SeqRecord child class with equality handling and prediction methods"""
    def __init__(self, seqrecord: SeqRecord, originhash = None):
        super().__init__(seq=seqrecord.seq,
            id=seqrecord.id,
            name=seqrecord.name,
            description=seqrecord.description,
            dbxrefs = seqrecord.dbxrefs,
            features = seqrecord.features,
            annotations = seqrecord.annotations,
            letter_annotations = seqrecord.letter_annotations)

        self.formatted_name = None

        self._selected_predictions = None

        self._ssrs = None
        self._srss = None
        self._rmds = None

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

    def call_predictions(self, strategy, is_circular=False):
        seq = str(self.seq.strip("\n").upper().replace("U", "T"))
        ssr_df, srs_df, rmd_df = predict(seq, strategy, is_circular)
        self._ssrs, self._srss, self._rmds = post_process(ssr_df,
                                                             srs_df,
                                                             rmd_df,
                                                             self,
                                                             is_circular)

    def from_same_file(self, other):
        if not self._originhash or not other._filehash:
            return False
        elif self._originhash != other._filehash:
            return False
        else:
            return True
