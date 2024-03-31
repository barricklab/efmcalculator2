import logging
from .short_seq_scan import scan_short_sequence
from .mut_rate_finder import get_mut_rate,get_recombo_rate

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())

MIN_SHORT_SEQ_LEN = 5
MAX_SHORT_SEQ_LEN = 17
UNKNOWN_REC_TYPE = "unknown"
SUB_RATE = float(2.2 * 10 ** (-10))



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
        
        
        first_find = True    
        loop_end = False
        ssr_count = 0
        sav_seq_attr = []
        mu_rate = 0
        for seq_attr in scan_short_sequence(seq, sub_seq, seq_len, isCircular, count):
            logger.info("test sequence finder")
            if seq_attr.length >= 2 and seq_attr.length <= 3:
                 _find_ssr(df,seq_attr,sav_seq_attr,first_find,loop_end,ssr_count)

            elif (seq_attr.length >= 6 and seq_attr.length <= 14):
                 _write_to_df(df,seq_attr,"SRS",str(count),mu_rate)

            elif seq_attr.length >= 15 and get_recombo_rate(seq_attr.length, seq_attr.distance,"ecoli" )!= 0:
                  _write_to_df(df,seq_attr,"RMD",str(count),get_recombo_rate(seq_attr.length, seq_attr.distance,"ecoli" ))

            if count <= 1 or str(sub_seq) not in visited_sequences:
                visited_sequences.add(str(sub_seq))

        if seq_attr.length >= 2 and seq_attr.length <= 3:
            loop_end = True
            _find_ssr(df,seq_attr)
        



def _find_ssr(df,seq_attr,first_find,loop_end,ssr_count,sav_seq_attr):
       if first_find == True:
           sav_seq_attr = seq_attr
           first_find = False
       if loop_end == False and seq_attr.distance == 1:
            ssr_count += 1
            first_find = False
       else:
            if ssr_count > 1:
                ssr_count += 1
                _write_to_df(df,sav_seq_attr,
                            "SSR", 
                            str(ssr_count), 
                            get_mut_rate(ssr_count,sav_seq_attr.length, "ecoli"))
            if sav_seq_attr.sub_seq == seq_attr.sub_seq:
                sav_seq_attr = seq_attr
                ssr_count = 0


def _write_to_df(df,seq_attr,rep_type, rep_rate, mu_rate):

    logger.debug (
        """Sequence = {} 
        : Size = {} 
        : Distance = {} 
        : Start-Pos = {} 
        : End-Pos = {}
        : rep_type = {}
        : rep_rate  = {}
        : mu_rate""".format(
        seq_attr.sub_seq,
        seq_attr.length,
        seq_attr.distance,
        seq_attr.start_pos,
        seq_attr.end_pos,
        rep_type,
        rep_rate,
        mu_rate
        ) 
    )
    df.loc[seq_attr.sub_seq + "~" + str(seq_attr.start_pos)] = [
                seq_attr.sub_seq,
                seq_attr.length,
                seq_attr.distance,
                seq_attr.start_pos,
                seq_attr.end_pos,
                rep_type,
                rep_rate,
                mu_rate
            ]


# EOF
