import polars as pl


def eval_top(ssr_df=None, srs_df=None, rmd_df=None, num_report: int = 10):
    # Get the top n=num_report contributors across all categories
    valid_dataframes = []
    if isinstance(ssr_df, pl.DataFrame) and not ssr_df.is_empty():
        ssr_df = ssr_df.sort(by="mutation_rate")
        ssr_df_top = (
            ssr_df.with_columns(source=pl.lit("SSR"))
            .sort(by="mutation_rate", descending=True)
            .select(pl.col(["repeat", "source", "mutation_rate", "predid"]))
            .head(num_report)
        )
        valid_dataframes.append(ssr_df_top)

    if isinstance(srs_df, pl.DataFrame) and not srs_df.is_empty():
        srs_df = srs_df.sort(by="mutation_rate")
        srs_df_top = (
            srs_df.with_columns(source=pl.lit("SRS"))
            .sort(by="mutation_rate", descending=True)
            .select(pl.col(["repeat", "source", "mutation_rate", "predid"]))
            .head(num_report)
        )
        valid_dataframes.append(srs_df_top)

    if isinstance(srs_df, pl.DataFrame) and not ssr_df.is_empty():
        rmd_df = rmd_df.sort(by="mutation_rate")
        rmd_df_top = (
            rmd_df.with_columns(source=pl.lit("RMD"))
            .sort(by="mutation_rate", descending=True)
            .select(pl.col(["repeat", "source", "mutation_rate", "predid"]))
            .head(num_report)
        )
        valid_dataframes.append(rmd_df_top)

    if len(valid_dataframes) > 0:
        merged_df = pl.concat(valid_dataframes)
        merged_df = merged_df.sort(by="mutation_rate", descending=True).head(num_report)
    # If no SSR, SRS, or RMD are found
    else:
        merged_df = srs_df

    return merged_df
