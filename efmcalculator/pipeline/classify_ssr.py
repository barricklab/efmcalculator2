import polars as pl

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

    collapsed_ssrs = (
        collapsed_ssrs
        .with_columns(
            pl.col("differences").cast(pl.List(pl.Float64))
            )
        .with_columns(
            truth_table=(
                pl.col("differences").list.eval(
                    (pl.element() != 0).or_(pl.element().is_null())
                )
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
    if collapsed_ssrs.height > 0:
        # Correct for circular SSR
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
                        pl.col("count").list.slice(1, pl.col("count").list.len() - 2),
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
