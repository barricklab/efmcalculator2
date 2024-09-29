import logging
import polars as pl
import tempfile
import itertools
from progress.spinner import Spinner
from progress.bar import IncrementalBar
import numpy as np
from .constants import MAX_SHORT_SEQ_LEN, MIN_SSR_LEN, MIN_SHORT_SEQ_LEN


class SeqAttr:
    def __init__(self, sub_seq, distance, start_pos, end_pos, note):
        self.sub_seq = str(sub_seq)
        self.length = len(sub_seq)
        self.distance = distance
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.note = note


class FakeBar:
    @staticmethod
    def next():
        return

    @staticmethod
    def finish():
        return


logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())


def _pairwise_slips(polars_df, column, is_circular) -> pl.DataFrame:
    """Recieve a polars dataframe with column of [List[type]]
    Returns the dataframe back with a pairwise list of positions"
    
    input:
        polars_df (dataframe): polars dataframe containing all repeats with lengths less than 16 base pairs
        column (str): name of column of [List[type]] in dataframe that contains the positions of each repeat
        is_circular (boolean): whether the scanned sequence is circular or not

    returns:
        polars dataframe containing all repeats with length less than 16 base pairs with a pairwise list of positions
    """

    # Shamelessly adapted from https://stackoverflow.com/questions/77354114/how-would-i-generate-combinations-of-items-within-polars-using-the-native-expres


    def map_function(list_o_things):
        return [
            sorted((thing_1, thing_2))
            for thing_1, thing_2 in itertools.combinations(list_o_things, 2)
        ]

    # Use linear pairing to find SSRs below the minimum SRS length (much faster)
    polars_df = polars_df.with_columns(
        pl.col('repeat').str.len_chars().alias('length')
    )

    linear_subset = polars_df.filter(pl.col('length') < MIN_SSR_LEN)
    pairwise = polars_df.filter(pl.col('length') >= MIN_SSR_LEN)

    linear_subset = _linear_slips(linear_subset, column, is_circular=is_circular)

    # Run pairwise for SRS's above minimum length
    pairwise = pairwise.lazy().with_columns(
        pl.col(column)
        .map_elements(map_function, return_dtype=pl.List(pl.List(pl.Int64)))
        .alias(f"pairings")
    )

    linear_subset = linear_subset
    pairwise = pairwise.collect().select(pl.col('repeat'), pl.col('pairings'))

    pairwise = pl.concat([linear_subset, pairwise])

    return pairwise


def _linear_slips(polars_df, column, is_circular=False) -> pl.LazyFrame:
    """Recieve a polars dataframe with column of [List[type]]
    Returns the dataframe back with a linear list of pairings
    
    input:
        polars_df (dataframe): polars dataframe containing all repeats with lengths less than 16 base pairs
        column (str): name of column of [List[type]] in dataframe that contains the positions of each repeat
        is_circular (boolean): whether the scanned sequence is circular or not

    returns:
        polars dataframe containing all repeats with length less than 16 base pairs with a linear list of pairings
    """

    nrows = polars_df.select(pl.len()).item()

    linear_df = (
        polars_df.with_columns(instances=pl.col(column).list.len())
        .with_columns(
            duplicate_column=pl.col(column)
            .list.tail(pl.col(column).list.len() - 1)
            .list.concat(pl.col(column).list.first())
        )
        .explode([column, "duplicate_column"])
        .select(
            [
                "repeat",
                pl.struct([column, "duplicate_column"]).alias("pairings"),
            ]
        )
        .group_by("repeat")
        .agg(pl.col("pairings"))
    )

    if not is_circular:
        linear_df = linear_df.with_columns(
            pairings=pl.col("pairings").list.head(pl.col("pairings").list.len() - 1)
        )

    linear_df = (
        linear_df.explode("pairings")
        .unnest("pairings")
        .select(
            pl.col("repeat"),
            pl.concat_list(pl.col("position"), pl.col("duplicate_column"))
            .list.sort()
            .alias("pairings"),
        )
        .group_by("repeat")
        .agg(pl.col("pairings").unique())
    )

    return linear_df


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
        distance=pl.col("pairings").list.diff().list.get(1) - pl.col("repeat_len")
    )

    if circular:
        distance_df = distance_df.with_columns(
            distance=pl.when(pl.col("distance") > seq_len / 2)
            .then(seq_len - pl.col("distance") - (2 * pl.col("repeat_len")))
            .otherwise(pl.col("distance")),
            wraparound=pl.when(pl.col("distance") > seq_len / 2)
            .then(True)
            .otherwise(False),
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


def _collapse_ssr(polars_df) -> pl.DataFrame:
    """Takes in dataframe of SSRs, returns dataframe of SSRs collapsed down.
    
    input:
        polars_df: polars datafrmae containing all repeats found, with a column for repeat type
    
    returns:
        polars dataframe containing only the SSRs found
    """

    collapsed_ssrs = (
        polars_df.filter(pl.col("category") == "SSR")
        .select(["repeat", "repeat_len", "pairings", "wraparound"])
        .lazy()
        # Collect the positions from all potentially participating SSRs
        .explode("pairings")
        .group_by(["repeat", "repeat_len"])
        .agg("pairings", "wraparound")
        .with_columns(
            positions=pl.col("pairings").list.unique().list.sort(),
        )
        .with_columns(differences=pl.col("positions").list.diff())
        .collect()
    )
    # Somehow couldnt figure out how to do this in pure polars

    # Identify start positions. If distance!=0, it a start of a repeat
    collapsed_ssrs = collapsed_ssrs.to_pandas()
    collapsed_ssrs["differences"] = (
        collapsed_ssrs["differences"] - collapsed_ssrs["repeat_len"]
    )
    collapsed_ssrs = pl.from_pandas(collapsed_ssrs)

    collapsed_ssrs = collapsed_ssrs.with_columns(
        truth_table=(
            pl.col("differences").list.eval(
                (pl.element() != 0).or_(pl.element().is_null())
            )
        )
    )

    # Apply the truth table
    collapsed_ssrs = collapsed_ssrs.to_pandas()
    collapsed_ssrs["starts"] = (collapsed_ssrs["truth_table"]) * (
        collapsed_ssrs["positions"] + 1
    )

    # Fill out
    collapsed_ssrs = pl.from_pandas(collapsed_ssrs)
    collapsed_ssrs = collapsed_ssrs.with_columns(
        starts=pl.col("starts")
        .explode()
        .replace(0, None)
        .forward_fill()
        .implode()
        .over("repeat")
    )

    # We had to apply an offeset with the truth table
    # To prevent bp-0 start positions from nulling out. This undoes that.
    collapsed_ssrs = collapsed_ssrs.to_pandas()
    collapsed_ssrs["starts"] = collapsed_ssrs["starts"] - 1
    collapsed_ssrs = pl.from_pandas(collapsed_ssrs).select(
        pl.col(["repeat", "repeat_len", "starts", "wraparound"])
    )

    # Count repeats from each start position
    collapsed_ssrs = (
        (
            collapsed_ssrs.lazy()
            .explode("starts")
            .group_by("repeat", "repeat_len", "starts", "wraparound")
            .count()
        )
        .rename({"starts": "start"})
        .collect()
    )

    # Correct for circular
    collapsed_ssrs = (
        collapsed_ssrs
        .sort("start")
        .group_by("repeat", "repeat_len")
        .agg(
            pl.col("start"),
            pl.col("count"), 
            pl.first("wraparound")
        )
        .with_columns(
            pl.when(pl.col("wraparound").list.contains(True))
            .then(
                pl.concat_list([
                    pl.col("count").list.slice(1, pl.col("count").list.lengths() - 2),
                    (pl.col("count").list.first() + pl.col("count").list.last())
                ]).alias("count")
            )
            .otherwise(pl.col("count")),
            pl.when(pl.col("wraparound").list.contains(True))
            .then(pl.col("start").list.slice(1).alias("start"))
            .otherwise(pl.col("start"))
        )
        .explode(
            pl.col("start"),
            pl.col("count")
        )
    )

    return collapsed_ssrs


def _build_seq_attr(sub_seq, seq_len, start_positions, isCircular, count):
    distance = 1
    start_pos = 0
    rem_start = 0
    end_pos = 0
    prv_end_pos = 0
    prv_start_pos = 0

    for start_pos in start_positions:
        note = ""
        end_pos = start_pos + len(sub_seq)
        # if end_pos >= seq_len:
        #   end_pos = end_pos - seq_len
        start_pos += 1
        distance = start_pos - prv_end_pos

        # fixes distance for
        if isCircular == True:
            if distance > seq_len / 2:
                distance = (seq_len + prv_start_pos) - end_pos
            if start_pos < seq_len and end_pos > seq_len:  # if repeat wraps around
                # fix end_pos
                end_pos = end_pos - seq_len
        if distance == 1:
            note = "SSR"
        # if overlapping
        if distance < 1:
            note = "skip for SSR"
            if count == 2:
                note = "skip"
        yield SeqAttr(sub_seq, distance, start_pos, end_pos, note)
        rem_start = start_pos + 1
        prv_end_pos = end_pos
        prv_start_pos = start_pos


def _find_repeat_positions(seq, sub_seq, seq_len, isCircular, count):
    """testing function for sequence scanning"""
    start_pos = 0
    rem_start = 0
    end_pos = 0
    prv_start_pos = 0
    start_positions = []

    start_pos = seq.find(sub_seq, rem_start)
    while start_pos != -1:
        start_positions.append(start_pos)
        start_pos += 1
        rem_start = start_pos + 1

        if len(sub_seq) == 1:
            start_pos = seq.find(sub_seq, rem_start - 1)
        else:
            start_pos = seq.find(sub_seq, rem_start)
    return start_positions
