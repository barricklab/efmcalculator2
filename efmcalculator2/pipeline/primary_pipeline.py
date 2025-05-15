import logging
import polars as pl
from collections import namedtuple
from collections import Counter, defaultdict
from rich import print
import Bio
import tempfile
from typing import List
from ..constants import MIN_SHORT_SEQ_LEN, MAX_SHORT_SEQ_LEN, UNKNOWN_REC_TYPE, SUB_RATE, MIN_SRS_LEN
from ..utilities import FakeBar

from .subsequence_curation import collect_subsequences, _scan_RMD, highly_mut
from .detection_strats import _pairwise_slips, _linear_slips
from .classify_ssr import _collapse_ssr
from .classify_srs_rmd import _calculate_distances, _categorize_efm


logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())


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

    # Remove SRS that are shorter than min_ssr_length
    repeat_df.filter(pl.col("repeat_len") >= MIN_SRS_LEN)

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
