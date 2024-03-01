import logging
from .short_seq_scan import scan_short_sequence

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())

MIN_SHORT_SEQ_LEN = 5
MAX_SHORT_SEQ_LEN = 17


def predict_RMDs(seq, df, seq_len, isCircular):
    _build_sub_seq_from_seq(seq, df, seq_len, isCircular)

def _build_sub_seq_from_seq(seq, df, seq_len, isCircular):
    for i, letter in enumerate(seq):
        for j in range(MAX_SHORT_SEQ_LEN, MIN_SHORT_SEQ_LEN, -1):
            if len(seq[i : i + j]) > MIN_SHORT_SEQ_LEN:
                sub_seq = seq[i : i + j]
                _find_short_seq(seq, sub_seq, df, seq_len, isCircular)


"""
    #circular short sequence look up
    k = MIN_SHORT_SEQ_LEN
    for j in (MIN_SHORT_SEQ_LEN, MAX_SHORT_SEQ_LEN,1):
        if (j > k):
            for l in range(len(seq) % k):
                roll_arr = numpy.roll([c for c in seq],l)
                roll_seq = ''.join(str(x) for x in roll_arr)
                sub_seq = Seq(roll_seq[0:k])
                find_short_seq(seq,sub_seq,df)
"""


def _find_short_seq(seq, sub_seq, df, seq_len, isCircular):
    count = seq.count_overlap(sub_seq)
    visited_sequences = set()
    if count > 1 and str(sub_seq) not in visited_sequences:
        df["Sequence"] = df["Sequence"].astype(str)
        # Filter for largest SRS only
        # df_filter = df[df['Sequence'].str.contains(str(sub_seq))]
        # if df_filter.empty:

        for seq_attr in scan_short_sequence(
            seq, sub_seq, seq_len, isCircular, count
        ):
            logger.debug(
                "Sequence = {} : Size = {} : Distance = {} : Start-Pos = {} : End-Pos = {}".format(
                    seq_attr.sub_seq,
                    seq_attr.length,
                    seq_attr.distance,
                    seq_attr.start_pos,
                    seq_attr.end_pos,
                )
            )
            df.loc[seq_attr.sub_seq + "~" + str(seq_attr.start_pos)] = [
                seq_attr.sub_seq,
                seq_attr.length,
                seq_attr.distance,
                seq_attr.start_pos,
                seq_attr.end_pos,
            ]
            if count <= 1 or str(sub_seq) not in visited_sequences:
                visited_sequences.add(str(sub_seq))

# EOF
