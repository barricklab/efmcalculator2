import logging
import polars as pl
from .short_seq_scan import _find_repeat_positions, _build_seq_attr
from .mut_rate_finder import get_mut_rate,get_recombo_rate
from collections import namedtuple
from progress.bar import IncrementalBar
from collections import Counter, defaultdict
from rich import print
import Bio
import tempfile



logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())

MIN_SHORT_SEQ_LEN = 0
MAX_SHORT_SEQ_LEN = 15
UNKNOWN_REC_TYPE = "unknown"
SUB_RATE = float(2.2 * 10 ** (-10))

class FakeBar():
    def next():
        return
    def finish():
        return

def predict_RMDs(seq, df, seq_len, isCircular, threads):
    _build_sub_seq_from_seq(seq, df, seq_len, isCircular, threads)
    #df.sort_values(['Size'], ascending=[True]).to_csv(output, index=False)
    # get RIP
    # Assuming df is your DataFrame

    tot_ssr_mut_rate = df["Mutation Rate"].where(df["Classifier"]=='SSR').sum()
    tot_rmd_mut_rate = df["Mutation Rate"].where(df["Classifier"]=='RMD').sum()

    results = df["Classifier"].value_counts()

    result = _find_rip(tot_ssr_mut_rate , tot_rmd_mut_rate)
    logger.info(f"RIP Score: {result['rip']} \n ------------------ \nRMDs: {result['rmd_sum']} \nSSRs: {result['ssr_sum']}, \nBase: {result['bps_sum']}")
    print(results)


def collect_subsequences(seq, window_max=16) -> pl.LazyFrame:
    '''Scans across a given input sequence and returns a list of subsequences'''
    if logger.isEnabledFor(logging.INFO):
        bar = IncrementalBar('Scanning for repeats', max=len(seq))
    else:
        bar = FakeBar()


    def scan_genome():
        # Probably room for optimizations here
        for i, _ in enumerate(seq):
            for j in range(MIN_SHORT_SEQ_LEN, MAX_SHORT_SEQ_LEN):
                if len(seq[i : i + j]) > MIN_SHORT_SEQ_LEN:
                    sub_seq = seq[i : i + j]
                    yield {'repeat': str(sub_seq), 'position': i}
            bar.next()

    repeats = (pl.LazyFrame(scan_genome()).
                groupby('repeat').
                agg(pl.col("position")).
                collect())
    bar.finish()

    return repeats



def _build_sub_seq_from_seq(seq, df, seq_len, isCircular, threads):

    # Curate target sequences

    repeat_df = collect_subsequences(seq)

    # Filter out sequences
    repeat_df = repeat_df.filter(pl.col('position').list.len()> 1)
    num_repeated_sequences = repeat_df.select(pl.len()).item()


    # --- Debugging
    def map_function(row):
        result = _find_repeat_positions(seq, row[0], seq_len, isCircular, len(row[1]))
        return (row[0], result)

    repeat_df = repeat_df.map_rows(map_function)
    repeat_df = repeat_df.rename({"column_0": 'repeat',
                    "column_1": 'position_corrected'})
    debug = repeat_df.join(repeat_df, on='repeat')
    for row in debug.rows(named=True):
        try:
            assert row['position'] == row['position_corrected']
        except:
            #print(row)
            pass

    repeat_df = debug
    position_source = 2 # 1 for rhobust, 2 for preexisting

    # -- end debugging


    # Get length of each repeat
    repeat_df = repeat_df.with_columns(pl.col("repeat").str.len_chars().alias("repeat_len"))
    #print(repeat_df.head())

    if logger.isEnabledFor(logging.INFO):
        bar = IncrementalBar('Calculating mutation rates', max=num_repeated_sequences)
    else:
        bar = FakeBar()

    def map_function(row):
        _find_short_seq(seq, row[0], df, seq_len, row[position_source], isCircular, bar)
        return ()
    repeat_df.map_rows(map_function)
    bar.finish()



def _find_short_seq(seq, sub_seq, df, seq_len, start_positions, isCircular, bar):
    sub_seq = Bio.Seq.Seq(data=sub_seq)
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

        for seq_attr in _build_seq_attr(sub_seq, seq_len, start_positions, isCircular, count):
            if seq_attr.length <= 5:
                ssrStruct= _find_ssr(df,seq_attr,ssrStruct.first_find,ssrStruct.loop_end,ssrStruct.ssr_count,ssrStruct.sav_seq_attr)

            elif (seq_attr.length >= 6 and seq_attr.length <= 15) and seq_attr.note != "skip":
                 if seq_attr.note == "SSR":
                      ssrStruct= _find_ssr(df,seq_attr,ssrStruct.first_find,ssrStruct.loop_end,ssrStruct.ssr_count,ssrStruct.sav_seq_attr)
                 else:
                      _write_to_df(df,seq_attr,"SRS",str(count),mu_rate)

            elif seq_attr.length >= 16 and seq_attr.note != "skip":
                  _write_to_df(df,seq_attr,"RMD",str(count),get_recombo_rate(seq_attr.length, seq_attr.distance,"ecoli" ))


        if seq_attr.length <= 5:
            loop_end = True
            _find_ssr(df,seq_attr,ssrStruct.first_find,loop_end,ssrStruct.ssr_count,ssrStruct.sav_seq_attr)

        elif (seq_attr.length >= 6 and seq_attr.length <= 15) and seq_attr.note == "SSR":
             loop_end = True
             _find_ssr(df,seq_attr,ssrStruct.first_find,loop_end,ssrStruct.ssr_count+1,ssrStruct.sav_seq_attr)
    bar.next()

        


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
                mut_rate = get_mut_rate(ssr_count,sav_seq_attr.length, "ecoli")
                _write_to_df(df,sav_seq_attr,
                            "SSR", 
                            str(ssr_count), 
                            mut_rate)
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
    #print(ssr_sum)
    #print(rmd_sum)
    base_rate = float(1000) * float(SUB_RATE)
    # Add in the mutation rate of an individual nucleotide
    r_sum = ssr_sum + rmd_sum + base_rate
    # Set the maximum rate sum to 1 for now.
    if r_sum > 1:
        r_sum = float(1)
    rel_rate = (float(r_sum) / float(base_rate))

    return {'rip': rel_rate, 'ssr_sum': ssr_sum, 'rmd_sum': rmd_sum, 'bps_sum': base_rate}
