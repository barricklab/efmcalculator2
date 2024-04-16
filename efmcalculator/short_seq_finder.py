import logging
from .short_seq_scan import scan_short_sequence
from .mut_rate_finder import get_mut_rate,get_recombo_rate
from collections import namedtuple

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())

MIN_SHORT_SEQ_LEN = 0
MAX_SHORT_SEQ_LEN = 17
UNKNOWN_REC_TYPE = "unknown"
SUB_RATE = float(2.2 * 10 ** (-10))



def predict_RMDs(seq, df, seq_len, isCircular):
    _build_sub_seq_from_seq(seq, df, seq_len, isCircular)
    #df.sort_values(['Size'], ascending=[True]).to_csv(output, index=False)
    # get RIP
    # Assuming df is your DataFrame
    tot_ssr_mut_rate = df["Mutation Rate"].where(df["Classifier"]=='SSR').sum()
    tot_rmd_mut_rate = df["Mutation Rate"].where(df["Classifier"]=='RMD').sum()

    logger.info(_find_rip(tot_ssr_mut_rate , tot_rmd_mut_rate ))

def _build_sub_seq_from_seq(seq, df, seq_len, isCircular):
    for i, letter in enumerate(seq):
        for j in range(MAX_SHORT_SEQ_LEN, MIN_SHORT_SEQ_LEN, -1):
            if len(seq[i : i + j]) > MIN_SHORT_SEQ_LEN:
                sub_seq = seq[i : i + j]
                _find_short_seq(seq, sub_seq, df, seq_len, isCircular)



def _find_short_seq(seq, sub_seq, df, seq_len, isCircular):
    count = seq.count_overlap(sub_seq)
    visited_sequences = set()
    if count > 1 and str(sub_seq) not in visited_sequences:
        df["Sequence"] = df["Sequence"].astype(str)
        
        
        ssrStruct = namedtuple("ssrStruct",["first_find", "loop_end", "ssr_count", "sav_seq_attr"])
        ssrStruct.first_find=True
        ssrStruct.loop_end=False
        ssrStruct.ssr_count=0
        ssrStruct.sav_seq_attr=[]
        mu_rate = 0
        for seq_attr in scan_short_sequence(seq, sub_seq, seq_len, isCircular, count):
            if seq_attr.length <= 5:
                ssrStruct= _find_ssr(df,seq_attr,ssrStruct.first_find,ssrStruct.loop_end,ssrStruct.ssr_count,ssrStruct.sav_seq_attr)

            elif (seq_attr.length >= 6 and seq_attr.length <= 15) and seq_attr.note != "skip":
                 if seq_attr.note == "SSR":
                      ssrStruct= _find_ssr(df,seq_attr,ssrStruct.first_find,ssrStruct.loop_end,ssrStruct.ssr_count,ssrStruct.sav_seq_attr)
                 else:
                      _write_to_df(df,seq_attr,"SRS",str(count),mu_rate)

            elif seq_attr.length >= 16 and get_recombo_rate(seq_attr.length, seq_attr.distance,"ecoli" )!= 0 and seq_attr.note != "skip":
                  _write_to_df(df,seq_attr,"RMD",str(count),get_recombo_rate(seq_attr.length, seq_attr.distance,"ecoli" ))

            if count <= 1 or str(sub_seq) not in visited_sequences:
                visited_sequences.add(str(sub_seq))

        if seq_attr.length <= 5:
            loop_end = True
            _find_ssr(df,seq_attr,ssrStruct.first_find,loop_end,ssrStruct.ssr_count,ssrStruct.sav_seq_attr)
        elif (seq_attr.length >= 6 and seq_attr.length <= 15) and seq_attr.note == "SSR":
             loop_end = True
             _find_ssr(df,seq_attr,ssrStruct.first_find,loop_end,ssrStruct.ssr_count+1,ssrStruct.sav_seq_attr)

        


def _find_ssr(df,seq_attr,first_find,loop_end,ssr_count,sav_seq_attr):
    if seq_attr.note != ("skip for SSR" or "skip"):
        if first_find == True:
            sav_seq_attr = seq_attr
            first_find = False
        if loop_end == False and seq_attr.distance == 1:
            ssr_count += 1
            first_find = False
        else:
            if (seq_attr.length >= 2 and ssr_count >= 3) or (seq_attr.length == 1 and ssr_count >= 4):
                _write_to_df(df,sav_seq_attr,
                            "SSR", 
                            str(ssr_count), 
                            get_mut_rate(ssr_count,sav_seq_attr.length, "ecoli"))
            if sav_seq_attr.sub_seq == seq_attr.sub_seq:
                sav_seq_attr = seq_attr
                ssr_count = 1
        
    ssrStruct = namedtuple("ssrStruct",["first_find", "loop_end", "ssr_count", "sav_seq_attr"])
    return ssrStruct(first_find,loop_end,ssr_count,sav_seq_attr)



def _write_to_df(df,seq_attr,rep_type, rep_rate, mu_rate):

    logger.debug ("Sequence = {}  : Size = {}  : Distance = {}  : Start-Pos = {} : End-Pos = {}: rep_type = {}: rep_rate  = {} : mu_rate".format(
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


def _find_rip(ssr_sum, rmd_sum):
    """
    Calculates an RIP score for given sequence
    :param repeats: List of dictionaries of repeats
    :param seq_len: Length of input sequence
    :return: Total predicted RIP score for whole sequence
    """ 
    print(ssr_sum)
    print(rmd_sum)
    base_rate = float(1000) * float(SUB_RATE)
    # Add in the mutation rate of an individual nucleotide
    r_sum = ssr_sum + rmd_sum + base_rate
    # Set the maximum rate sum to 1 for now.
    if r_sum > 1:
        r_sum = float(1)
    rel_rate = (float(r_sum) / float(base_rate))

    return {'rip': rel_rate, 'ssr_sum': ssr_sum, 'rmd_sum': rmd_sum, 'bps_sum': base_rate}
