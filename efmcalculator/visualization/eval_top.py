import polars as pl


def eval_top(ssr_df=None, srs_df=None, rmd_df=None, num_report: int = 5):
    # Get the top n=num_report contributors across all categories
    valid_dataframes = []
    if isinstance(ssr_df, pl.DataFrame) and not ssr_df.is_empty():
        ssr_df = ssr_df.sort(by="mutation_rate")
        ssr_df_top = (
            ssr_df.with_columns(source=pl.lit("SSR"))
            .sort(by="mutation_rate", descending=True)
            .select(pl.col(["source", "mutation_rate"]))
            .head(num_report)
        )
        valid_dataframes.append(ssr_df_top)

    if isinstance(srs_df, pl.DataFrame) and not srs_df.is_empty():
        srs_df = srs_df.sort(by="mutation_rate")
        srs_df_top = (
            srs_df.with_columns(source=pl.lit("SRS"))
            .sort(by="mutation_rate", descending=True)
            .select(pl.col(["source", "mutation_rate"]))
            .head(num_report)
        )
        valid_dataframes.append(srs_df_top)

    if isinstance(srs_df, pl.DataFrame) and not ssr_df.is_empty():
        rmd_df = rmd_df.sort(by="mutation_rate")
        rmd_df_top = (
            rmd_df.with_columns(source=pl.lit("RMD"))
            .sort(by="mutation_rate", descending=True)
            .select(pl.col(["source", "mutation_rate"]))
            .head(num_report)
        )
        valid_dataframes.append(rmd_df_top)

    merged_df = pl.concat(valid_dataframes)
    merged_df = merged_df.sort(by="mutation_rate", descending=True).head(num_report)
    count = merged_df.select(pl.col("source").value_counts()).unnest("source").to_dict()
    count = dict(zip(count["source"], count["count"]))

    # Apply this information to the dataframes

    if isinstance(ssr_df, pl.DataFrame) and not ssr_df.is_empty():
        ssr_df = (
            ssr_df.sort(by="mutation_rate", descending=True)
            .with_row_index()
            .with_columns(
                show=pl.when(pl.col("index") < count.get("SSR", 0))
                .then(pl.lit(True))
                .otherwise(False)
            )
        )
    if isinstance(srs_df, pl.DataFrame) and not srs_df.is_empty():
        srs_df = (
            srs_df.sort(by="mutation_rate", descending=True)
            .with_row_index()
            .with_columns(
                show=pl.when(pl.col("index") < count.get("SRS", 0))
                .then(pl.lit(True))
                .otherwise(False)
            )
        )
    if isinstance(rmd_df, pl.DataFrame) and not rmd_df.is_empty():
        rmd_df = (
            rmd_df.sort(by="mutation_rate", descending=True)
            .with_row_index()
            .with_columns(
                show=pl.when(pl.col("index") < count.get("RMD", 0))
                .then(pl.lit(True))
                .otherwise(False)
            )
        )

    return ssr_df, srs_df, rmd_df
