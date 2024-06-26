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
from typing import Set
from .constants import MIN_SHORT_SEQ_LEN, MAX_SHORT_SEQ_LEN, UNKNOWN_REC_TYPE, SUB_RATE
from .utilities import FakeBar

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())


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
                    yield {"repeat": str(sub_seq), "position": i}
            bar.next()

    repeats = (
        pl.LazyFrame(scan_genome()).groupby("repeat").agg(pl.col("position")).collect()
    )
    bar.finish()

    return repeats


def _scan_RMD(df: pl.DataFrame) -> pl.DataFrame:
    """Scans for RMDs"""
    knwon_long_repeats = df.filter(pl.col("repeat_len") == MAX_SHORT_SEQ_LEN)

    # @TODO KEVIN - RMD repeat logic here

    return df  # In the same format as df alongside the original data


def predict(seq: str, strategy: str, isCircular: bool) -> Set[pl.DataFrame]:
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
    # -- end debugging

    # Create list of positions
    if strategy == "pairwise":
        repeat_df = _pairwise_slips(repeat_df, "position")
    elif strategy == "linear":
        repeat_df = _linear_slips(repeat_df, "position", isCircular)
    else:
        raise ValueError("Invalid strategy")

    repeat_df = repeat_df.drop("position")
    repeat_df = repeat_df.drop("position_corrected")
    repeat_df = repeat_df.explode("pairings")

    # Get length of each repeat
    repeat_df = repeat_df.with_columns(
        pl.col("repeat").str.len_chars().alias("repeat_len")
    )

    # Upgrade long SRSs to RMDs
    repeat_df = _scan_RMD(repeat_df)

    # Calculate Distances
    repeat_df = _calculate_distances(repeat_df, seq_len, isCircular)
    repeat_df = repeat_df.filter(pl.col("distance") >= 0)

    # Categorize positions
    repeat_df = _categorize_efm(repeat_df).collect()
    print(repeat_df)

    # Collapse SSRs down
    ssr_df = _collapse_ssr(repeat_df).select(
        pl.col(["repeat", "repeat_len", "start", "count"])
    )

    srs_df = (
        repeat_df.filter(pl.col("category") == "SRS")
        .select(pl.col(["repeat", "repeat_len", "pairings", "distance"]))
        .rename({"pairings": "positions"})
    )

    rmd_df = repeat_df.filter(pl.col("category") == "RMD").select(
        pl.col(["repeat", "repeat_len", "pairings", "distance"])
    )
    rmd_df = rmd_df.rename({"pairings": "positions"})

    return (ssr_df, srs_df, rmd_df)
