import polars as pl


def filter_ssrs(ssr_dataframe):
    ssr_dataframe = (
        ssr_dataframe.lazy()
        # Filter based on SSR definition
        .filter(
            (pl.col("repeat_len") >= 2)
            .and_(pl.col("count") >= 3)
            .or_((pl.col("repeat_len") == 1).and_(pl.col("count") >= 4))
        )
    ).collect()

    # Keep SSR with lowest length
    ssr_dataframe = (
        ssr_dataframe.sort("repeat_len")
        .group_by(pl.col("start"))
        .head(1)
    )

    # Delete SSR within other more important SSRs ("CACACA" will have SSR of "CA" and "AC")
    ssr_dataframe = (
    ssr_dataframe.sort(pl.col("start"))
    .with_columns(
        pl.col("repeat").shift(1).str.tail(1).alias("last_tail"),
        pl.col("repeat").str.head(1).alias("head"),
        pl.col("start").shift(1).alias("last_start"),
        )
    .filter(
        (pl.col("last_start").is_null()) |
        (pl.col("start") != pl.col("last_start") + 1) |
        (pl.col("last_tail") != pl.col("head"))
        )
    .select(["repeat", "repeat_len", "start", "count"])
    )

    return ssr_dataframe


def filter_rmds(rmd_dataframe):
    # Delete redundant SRS repeats @TODO

    return rmd_dataframe


# No longer necessary filters (already covered by pairwise approach)
# - Delete rows with overlapping repeats and only 2 occurrences (not a real repeat)
