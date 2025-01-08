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
from collections import namedtuple
from progress.bar import IncrementalBar
from collections import Counter, defaultdict
from rich import print
import Bio
import tempfile
from typing import List
from .constants import MIN_SHORT_SEQ_LEN, MAX_SHORT_SEQ_LEN, UNKNOWN_REC_TYPE, SUB_RATE
from .utilities import FakeBar
from .bad_state_mitigation import detect_special_cases
import streamlit as st

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())


def collect_subsequences(seq, isCircular, window_max=16) -> pl.LazyFrame:
    """Scans across a given input sequence and returns a list of subsequences

    input:
        seq (string): the sequence to be scanned
        isCircular (boolean): Whether the sequence is circular or not

    returns:
        polars dataframe containing every repeat smaller than 16 base pairs in the input sequence

    """

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
                    yield {"repeat": str(sub_seq), "position": i}
            bar.next()

    repeats = (
        pl.LazyFrame(scan_genome()).group_by("repeat").agg(pl.col("position")).collect()
    ).cast({"position": pl.List(pl.Int32)})

    bar.finish()

    return repeats


def _scan_RMD(df: pl.DataFrame, seq, seq_len, isCircular) -> pl.DataFrame:
    """Scans for RMDs

    input:
        df (dataframe): dataframe containing all repeats smaller than 16 base pairs
        seq (string): the sequence to be scanned
        seq_len (int): the length of the sequence
        isCircular (boolean): whether the sequence is circular or not

    returns:
        polars dataframe containing every repeat in the input sequence

    """

    known_long_repeats = df.filter(pl.col("repeat_len") == (MAX_SHORT_SEQ_LEN - 1))
    RMD_df = pl.DataFrame(
        {
            "repeat": pl.Series("repeat", [], pl.Utf8),
            "pairings": pl.Series("pairings", [], pl.List(pl.Int64)),
            "repeat_len": pl.Series("repeat_len", [], pl.Int64),
        }
    )

    def check_larger_repeats(positions, seq):
        completed = False
        step = 50
        length = MAX_SHORT_SEQ_LEN
        pos1 = positions[0]
        pos2 = positions[1]
        largest = False
        wrap = 0

        while not completed:
            prvlength = length
            length += step
            # already know they are 15 bp repeats
            if prvlength < 16:
                prvlength = 16

            # pos2 is always after pos1
            if pos2 + length > seq_len:
                if isCircular:
                    wrap += step
                    seq = seq + seq[wrap - step : wrap]
                else:
                    length = seq_len - pos2
                    largest = True
            # largest possible repeat is seq_len/2 bp
            if length >= int(seq_len / 2):
                largest = True
                length = int(seq_len / 2)
            # prevent out of bound error
            if len(seq) - pos2 < length:
                length = len(seq) - pos2

            # if sequences not equal, then largest repeat has been passed
            if (seq[pos1 : pos1 + length] != seq[pos2 : pos2 + length]) or (largest):
                # iterate 1 by 1.
                # Uses length+1 because in range stops before last index
                for j in range(prvlength, length + 1):
                    # uses j-1 because substrings end before last index
                    if seq[pos1 + (j - 1)] == seq[pos2 + (j - 1)]:
                        sub_seq = seq[pos1 : pos1 + j]
                        yield {
                            "repeat": str(sub_seq),
                            "pairings": [pos1, pos2],
                            "repeat_len": j,
                        }
                    else:
                        break
                completed = True

    def store_RMD(positions, seq):
        nonlocal RMD_df
        repeats = pl.DataFrame(check_larger_repeats(positions, seq))

        # if larger repeats were found, then repeats will have 3 columns ("repeat", "pairings", "repeat_len")
        if repeats.width == 3:
            RMD_df = pl.concat([RMD_df, repeats])
        return None

    # Apply the function to the DataFrame
    known_long_repeats.with_columns(
        pl.col("pairings").map_elements(
            lambda pairings: store_RMD(pairings, seq), return_dtype=pl.List(pl.Null)
        )
    )
    RMD_df = RMD_df.with_columns(
        pl.col("pairings").cast(pl.List(pl.Int32)),
        pl.col("repeat_len").cast(pl.Int32))

    df = pl.concat([df, RMD_df])

    return df  # In the same format as df alongside the original data

def predict(seq: str, strategy: str, isCircular: bool) -> List[pl.DataFrame]:
    """Scans and predicts SSRs and RMDs. Returns dataframes representing each

    input:
        seq (string): the sequence to be scanned
        strategy (str): the scanning strategy to be used. Either pairwise or linear
        isCircular (boolean): whether the sequence is circular or not

    returns:
        ssr_df (dataframe): dataframe containing all the SSRs found in the sequence
        srs_df (dataframe): dataframe containing all the SRSs found in the sequence
        rmd_df (dataframe): dataframe containing all the RMDs found in the sequence

    """
    seq = seq.upper().replace(" ", "")
    seq_len = len(seq)

    valid_strategies = ["pairwise", "linear"]
    if strategy not in valid_strategies:
        raise ValueError(
            f"Invalid strategy: {strategy}. Must be one of {valid_strategies}"
        )

    detect_special_cases(seq, circular=isCircular)

    # Curate target sequences

    repeat_df = collect_subsequences(seq, isCircular)

    # Filter out sequences
    repeat_df = repeat_df.filter(pl.col("position").list.len() > 1)
    num_repeated_sequences = repeat_df.select(pl.len()).item()
    # -- end debugging

    # Create list of positions
    if strategy == "pairwise":
        repeat_df = _pairwise_slips(repeat_df, "position", isCircular)
    elif strategy == "linear":
        repeat_df = _linear_slips(repeat_df, "position", isCircular)
    else:
        raise ValueError("Invalid strategy")

    if "position" in repeat_df:
        repeat_df = repeat_df.drop("position")
    if "position_corrected" in repeat_df:
        repeat_df = repeat_df.drop("position_corrected")
    repeat_df = repeat_df.explode("pairings")
    # Get length of each repeat
    repeat_df = repeat_df.with_columns(
        pl.col("repeat").str.len_chars().alias("repeat_len").cast(pl.Int32)
    )


    # Upgrade long SRSs to RMDs
    repeat_df = _scan_RMD(repeat_df, seq, seq_len, isCircular)

    # Calculate Distances
    repeat_df = _calculate_distances(repeat_df, seq_len, isCircular)
    repeat_df = repeat_df.filter(pl.col("distance") >= 0)
    repeat_df = repeat_df.unique()

    # Categorize positions
    repeat_df = _categorize_efm(repeat_df)
    # Collapse SSRs down
    ssr_df = _collapse_ssr(repeat_df).select(
        pl.col(["repeat", "repeat_len", "start", "count"])
    )

    # Process and Split SRS and RMD

    repeat_df = repeat_df.lazy().filter(pl.col("category") != "SSR").collect()
    if len(repeat_df) > 0:
        repeat_df = (
            repeat_df.lazy()
            .select(
                pl.col(["repeat", "repeat_len", "pairings", "distance", "category"])
            )
            .collect()
            .rechunk()  # Weird issue with invalid pairing state
            .lazy()
            .with_columns(
                pl.col("pairings").list.to_struct(
                    fields=[
                        "first_repeat",
                        "second_repeat",
                    ]
                )
            )
            .unnest("pairings")
        ).collect()
    else:
        schema = {
            "repeat": pl.Utf8,
            "repeat_len": pl.Int32,
            "first_repeat": pl.Int32,
            "second_repeat": pl.Int32,
            "distance": pl.Int32,
            "category": pl.Categorical,
        }
        repeat_df = pl.DataFrame(schema=schema)

    srs_df = repeat_df.filter(pl.col("category") == "SRS").select(
        pl.col(["repeat", "repeat_len", "first_repeat", "second_repeat", "distance"])
    )
    rmd_df = repeat_df.filter(pl.col("category") == "RMD").select(
        pl.col(["repeat", "repeat_len", "first_repeat", "second_repeat", "distance"])
    )

    return [ssr_df, srs_df, rmd_df]
