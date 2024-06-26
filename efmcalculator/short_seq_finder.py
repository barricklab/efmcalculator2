import logging
import polars as pl
from .short_seq_scan import (
    _find_repeat_positions,
    _build_seq_attr,
    _linear_slips,
    _pairwise_slips,
    _calculate_distances,
    _categorize_efm,
    _collapse_ssr,
)
from .mut_rate_finder import get_mut_rate, get_recombo_rate
from collections import namedtuple
from progress.bar import IncrementalBar
from collections import Counter, defaultdict
from rich import print
import Bio
import tempfile


logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())

MIN_SHORT_SEQ_LEN = 0
MAX_SHORT_SEQ_LEN = 16
UNKNOWN_REC_TYPE = "unknown"
SUB_RATE = float(2.2 * 10 ** (-10))


class FakeBar:
    def next():
        return

    def finish():
        return


"""
def predict(seq, df, isCircular, threads):
    _build_sub_seq_from_seq(seq, df, isCircular, threads)
    #df.sort_values(['Size'], ascending=[True]).to_csv(output, index=False)
    # get RIP
    # Assuming df is your DataFrame

    tot_ssr_mut_rate = df["Mutation Rate"].where(df["Classifier"]=='SSR').sum()
    tot_rmd_mut_rate = df["Mutation Rate"].where(df["Classifier"]=='RMD').sum()

    results = df["Classifier"].value_counts()

    result = _find_rip(tot_ssr_mut_rate , tot_rmd_mut_rate)
    logger.info(f"RIP Score: {result['rip']} \n ------------------ \nRMDs: {result['rmd_sum']} \nSSRs: {result['ssr_sum']}, \nBase: {result['bps_sum']}")
    print(results)
"""


def collect_subsequences(seq, isCircular, window_max=16) -> pl.LazyFrame:
    """Scans across a given input sequence and returns a list of subsequences"""
    if logger.isEnabledFor(logging.INFO):
        bar = IncrementalBar("Scanning for repeats", max=len(seq))
    else:
        bar = FakeBar()

    seq_len = len(seq)

    if isCircular:
        # adds first 20 bp to end
        seq = seq + seq[0:20]

    def scan_genome():
        # Probably room for optimizations here
        for i, _ in enumerate(seq[:seq_len]):
            for j in range(MIN_SHORT_SEQ_LEN, MAX_SHORT_SEQ_LEN):
                if len(seq[i : i + j]) > MIN_SHORT_SEQ_LEN:
                    sub_seq = seq[i : i + j]
                    if seq.count_overlap(sub_seq) > 1: 
                        yield {'repeat': str(sub_seq), 'position': i}
                        repeats_left = seq.count_overlap(sub_seq)
                    else:
                        repeats_left = 0
                        break
            
            # scan from smallest RMD to largest possible RMD, iterating by 50
            step = 50
            j = MAX_SHORT_SEQ_LEN
            max_repeat_len = int(seq_len/2)

            while j <= max_repeat_len:
                # if repeats_left <= 1, then all repeats have been found already 
                if repeats_left > 1:
                    sub_seq = seq[i : i  + j]
                    # seq.count_overlap < repeats left means largest repeat is between previous length and current length
                    # if j == max_repeat_len, then the largest repeat must be between previous length and current length
                    if (seq.count_overlap(sub_seq) < repeats_left) or (j == max_repeat_len):
                        start = j - step
                        # repeats smaller than MAX_SHORT_SEQ_LEN have already been found
                        if start < MAX_SHORT_SEQ_LEN:
                            start = MAX_SHORT_SEQ_LEN

                        # start iterating by 1 from previous length to find largest repeat
                        for k in range(start, j + 1):
                            sub_seq = seq[i : i + k]
                            if seq.count_overlap(sub_seq) > 1:
                                yield {'repeat': str(sub_seq), 'position': i}
                            # no more repeats
                            else:
                                break
                        # find how many repeats are still left
                        repeats_left = seq.count_overlap(sub_seq)
                    # changes j so that the largest repeat possible (seq_len/2) is always included
                    # if j == max_repeat_len, then this code has already ran and shouldn't be ran again
                    if j + step >= max_repeat_len and j != max_repeat_len:
                        j = max_repeat_len - step
                else:
                    break
                j += step

            bar.next()

    repeats = (
        pl.LazyFrame(scan_genome()).groupby("repeat").agg(pl.col("position")).collect()
    )
    bar.finish()

    return repeats


def predict(seq, df, strategy, isCircular, threads) -> (pl.DataFrame, pl.DataFrame):
    """Scans and predicts SSRs and RMDs. Returns dataframes representing each"""
    seq_len = len(seq)

    valid_strategies = ["pairwise", "linear"]
    if strategy not in valid_strategies:
        raise ValueError(
            f"Invalid strategy: {strategy}. Must be one of {valid_strategies}"
        )

    # Curate target sequences

    repeat_df = collect_subsequences(seq, isCircular)

    # Filter out sequences
    repeat_df = repeat_df.filter(pl.col("position").list.len() > 1)
    num_repeated_sequences = repeat_df.select(pl.len()).item()

    # --- Debugging
    """
    def map_function(row):
        result = _find_repeat_positions(seq, row[0], seq_len, isCircular, len(row[1]))
        return (row[0], result)

    debug = repeat_df.map_rows(map_function)
    debug = debug.rename({"column_0": 'repeat',
                    "column_1": 'position'})
    debug = debug.join(debug, on='repeat')
    debug = debug.rename({"position_right": 'position_corrected'})

    for row in debug.rows(named=True):
        try:
            assert row['position'] == row['position_corrected']
        except:
            print(row)
            pass

    repeat_df = debug
    #print(repeat_df)
    """
    position_source = 1  # 1 for rhobust, 2 for preexisting

    # -- end debugging

    # Create list of positions
    if strategy == "pairwise":
        repeat_df = _pairwise_slips(repeat_df, "position")
    elif strategy == "linear":
        repeat_df = _linear_slips(repeat_df, "position", isCircular)
    else:
        return

    repeat_df = repeat_df.drop("position")
    repeat_df = repeat_df.drop("position_corrected")
    repeat_df = repeat_df.explode("pairings")

    # Get length of each repeat
    repeat_df = repeat_df.with_columns(
        pl.col("repeat").str.len_chars().alias("repeat_len")
    )

    # Calculate Distances
    repeat_df = _calculate_distances(repeat_df, seq_len, isCircular)
    repeat_df = repeat_df.filter(pl.col("distance") >= 0)

    # Categorize positions
    repeat_df = _categorize_efm(repeat_df)

    # Collapse SSRs down
    ssr_df = _collapse_ssr(repeat_df).select(
        pl.col(["repeat", "repeat_len", "start", "count"])
    )

    rmd_df = repeat_df.filter(pl.col("is_RMD") == True).select(
        pl.col(["repeat", "repeat_len", "pairings", "distance"])
    )
    rmd_df = rmd_df.rename({"pairings": "positions"})

    # ---------------------- LEGACY CODE TO BE REMOVED

    # print(pairwise_df.filter(pl.col('is_RMD') == False))

    if logger.isEnabledFor(logging.INFO):
        bar = IncrementalBar("Calculating mutation rates", max=num_repeated_sequences)
    else:
        bar = FakeBar()

    def map_function(row):
        _find_short_seq(seq, row[0], df, seq_len, row[position_source], isCircular)
        bar.next()
        return ()

    # repeat_df.map_rows(map_function)
    bar.finish()

    # ---------------------- LEGACY CODE TO BE REMOVED

    return (ssr_df, rmd_df)


def _find_short_seq(seq, sub_seq, df, seq_len, start_positions, isCircular):
    sub_seq = Bio.Seq.Seq(data=sub_seq)
    count = seq.count_overlap(sub_seq)
    visited_sequences = set()
    if count > 1 and str(sub_seq) not in visited_sequences:
        df["Sequence"] = df["Sequence"].astype(str)

        ssrStruct = namedtuple(
            "ssrStruct", ["first_find", "loop_end", "ssr_count", "sav_seq_attr"]
        )
        ssrStruct.first_find = True
        ssrStruct.loop_end = False
        ssrStruct.ssr_count = 0
        ssrStruct.sav_seq_attr = []
        mu_rate = 0

        for seq_attr in _build_seq_attr(
            sub_seq, seq_len, start_positions, isCircular, count
        ):
            if seq_attr.length <= 5:
                ssrStruct = _find_ssr(
                    df,
                    seq_attr,
                    ssrStruct.first_find,
                    ssrStruct.loop_end,
                    ssrStruct.ssr_count,
                    ssrStruct.sav_seq_attr,
                )

            elif (
                seq_attr.length >= 6 and seq_attr.length <= 15
            ) and seq_attr.note != "skip":
                if seq_attr.note == "SSR":
                    ssrStruct = _find_ssr(
                        df,
                        seq_attr,
                        ssrStruct.first_find,
                        ssrStruct.loop_end,
                        ssrStruct.ssr_count,
                        ssrStruct.sav_seq_attr,
                    )
                else:
                    _write_to_df(df, seq_attr, "SRS", str(count), mu_rate)

            elif seq_attr.length >= 16 and seq_attr.note != "skip":
                _write_to_df(
                    df,
                    seq_attr,
                    "RMD",
                    str(count),
                    get_recombo_rate(seq_attr.length, seq_attr.distance, "ecoli"),
                )

        if seq_attr.length <= 5:
            loop_end = True
            _find_ssr(
                df,
                seq_attr,
                ssrStruct.first_find,
                loop_end,
                ssrStruct.ssr_count,
                ssrStruct.sav_seq_attr,
            )

        elif (
            seq_attr.length >= 6 and seq_attr.length <= 15
        ) and seq_attr.note == "SSR":
            loop_end = True
            _find_ssr(
                df,
                seq_attr,
                ssrStruct.first_find,
                loop_end,
                ssrStruct.ssr_count + 1,
                ssrStruct.sav_seq_attr,
            )


def _find_ssr(df, seq_attr, first_find, loop_end, ssr_count, sav_seq_attr):
    if seq_attr.note != ("skip for SSR" or "skip"):
        if first_find == True:
            sav_seq_attr = seq_attr
            first_find = False
        if loop_end == False and seq_attr.distance == 1:
            ssr_count += 1
            first_find = False
        else:
            if (seq_attr.length >= 2 and ssr_count >= 3) or (
                seq_attr.length == 1 and ssr_count >= 4
            ):
                mut_rate = get_mut_rate(ssr_count, sav_seq_attr.length, "ecoli")
                _write_to_df(df, sav_seq_attr, "SSR", str(ssr_count), mut_rate)
            if sav_seq_attr.sub_seq == seq_attr.sub_seq:
                sav_seq_attr = seq_attr
                ssr_count = 1

    ssrStruct = namedtuple(
        "ssrStruct", ["first_find", "loop_end", "ssr_count", "sav_seq_attr"]
    )
    return ssrStruct(first_find, loop_end, ssr_count, sav_seq_attr)


def _write_to_df(df, seq_attr, rep_type, rep_rate, mu_rate):
    logger.debug(
        "Sequence = {}  : Size = {}  : Distance = {}  : Start-Pos = {} : End-Pos = {}: rep_type = {}: rep_rate  = {} : mu_rate".format(
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
        mu_rate,
    ]
