import logging
import polars as pl
import tempfile
import itertools
from progress.spinner import Spinner
from progress.bar import IncrementalBar
import numpy as np
from ..constants import MAX_SHORT_SEQ_LEN, MIN_SHORT_SEQ_LEN


logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())


def _calculate_distances(polars_df, seq_len, circular) -> pl.LazyFrame:
    """
    calculates the distance between occurrences of each repeat

    input:
        polars_df: polars dataframe containing all repeats found in the scanned sequence
        seq_len (int): length of the scanned sequence
        circular: whether the scanned sequence is circular or not

    returns:
        inputted polars dataframe with a column of distances between each repeat
    """
    distance_df = polars_df.with_columns(
        distance=pl.col("pairings").list.diff().list.get(1) - pl.col("repeat_len").cast(pl.Int32)
    )
    if circular:
         distance_df = (
            distance_df
            .with_columns((pl.col("pairings").list.get(0) - (pl.col("pairings").list.get(1) + pl.col("repeat_len") - seq_len)).alias("alt_len"))
            .with_columns(
                distance=pl.when(pl.col("alt_len") >= 0)
                .then(
                    pl.when(pl.col("distance") > pl.col("alt_len"))
                    .then(pl.col("alt_len"))
                    .otherwise(pl.col("distance"))
                )
                .otherwise(pl.col("distance")).cast(pl.Int32),
                wraparound=pl.when(pl.col("alt_len") >= 0)
                .then(
                    pl.when(pl.col("distance") > pl.col("alt_len"))
                    .then(True)
                    .otherwise(False)
                )
                .otherwise(False)
                )
        )
    else:
        distance_df = distance_df.with_columns(
            wraparound = False,
        )
    return distance_df

def _categorize_efm(polars_df) -> pl.DataFrame:
    """
    categorizes every repeat as SSR, SRS, or RMD

    input:
        polars_df: polars dataframe containing all repeats found in the scanned sequence

    returns:
        inputted polars dataframe with a column for repeat type, either SSR, SRS, or RMD
    """
    categorized_df = polars_df.with_columns(
        category=pl.when(pl.col("distance") > 0)
        .then(
            pl.when(pl.col("repeat_len") > MIN_SHORT_SEQ_LEN)
            .then(
            pl.when(pl.col("repeat_len") < MAX_SHORT_SEQ_LEN)
            .then(pl.lit("SRS"))
            .otherwise(pl.lit("RMD"))
            )
        )
        .otherwise(pl.lit("SSR"))
        .cast(pl.Categorical)
    )
    return categorized_df
