import polars as pl
import itertools
from ..constants import MIN_SRS_LEN, MAX_SHORT_SEQ_LEN


def _pairwise_slips(polars_df, column, is_circular, copies_cap=40) -> pl.DataFrame:
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

    linear_subset = polars_df.filter(pl.col('length') < MIN_SRS_LEN)
    pairwise = polars_df.filter(pl.col('length') >= MIN_SRS_LEN)

    high_mut_df = (
        pairwise.filter(
            pl.col("length") == MAX_SHORT_SEQ_LEN-1,
            pl.col("position").list.len() >= copies_cap
            )
    )

    # terminate if there are 40+ occurrences of a single repeat
    if len(high_mut_df) > 0:
        raise ValueError("This sequence is highly mutagenic. Stopping execution")

    linear_subset = _linear_slips(linear_subset, column, is_circular=is_circular)

    # Run pairwise for SRS's above minimum length
    pairwise = pairwise.lazy().with_columns(
        pl.col(column)
        .map_elements(map_function, return_dtype=pl.List(pl.List(pl.Int64)))
        .alias(f"pairings").cast(pl.List(pl.List(pl.Int32)))
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
