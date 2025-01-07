import polars as pl


def filter_ssrs(ssr_dataframe, seq_len, circular):
    ssr_dataframe = (
        ssr_dataframe.lazy()
        # Filter based on SSR definition
        .filter(
            (pl.col("repeat_len") >= 3)
            .and_(pl.col("count") >= 3)
            .or_((pl.col("repeat_len") <= 2).and_(pl.col("count") >= 4))
        )   #! EFM1 is bugged with repeat_ len>= 3 and <= 2 bp
    ).collect()

    # Keep SSR with lowest length that start at same position
    ssr_dataframe = (
        ssr_dataframe.sort("repeat_len")
        .group_by(pl.col("start"))
        .head(1)
    )

    # Delete SSR within other more important SSRs ("CACACA" will have SSR of "CA" and "AC")
    ssr_dataframe = (
    ssr_dataframe.sort(pl.col("start"))
        .with_columns(
            pl.when(circular)
            .then(
                pl.when((pl.col("start") + (pl.col("count") * pl.col("repeat_len"))) >= seq_len)
                .then(pl.col("start")-seq_len)
                .otherwise(pl.col("start"))
                )
        .otherwise(pl.col("start"))
        )
        .with_columns(
        pl.col("start").shift(1).alias("last_start"),
        pl.col("count").shift(1).alias("last_count"),
        )

        .filter(
            # gets rid only if starts 1 bp after last SSR, and has lower or equal count
            (
                (pl.col("last_start").is_null()) |
                (pl.col("start") != pl.col("last_start") + 1) |
                (pl.col("count") > pl.col("last_count"))
            ) &
            # keep SSR with lowest length when not starting from same position
            (
                (pl.col("last_start").is_null()) |
                (pl.col("last_start") + pl.col("repeat_len") != pl.col("start"))
            )
        )
    .select(["repeat", "repeat_len", "start", "count"])
    # fix start values
    .with_columns(
        pl.when(pl.col("start") < 0)
        .then(pl.col("start") + seq_len)
        .otherwise(pl.col("start"))
        )
    )

    return ssr_dataframe


def filter_direct_repeats(rmd_dataframe, srs_dataframe, seq_len, ssr_dataframe, circular):
    # Delete redundant SRS repeats

    # label RMDs and SRSs and combine into 1 df
    rmd_dataframe = rmd_dataframe.with_columns(pl.lit("RMD").alias("type"))

    srs_dataframe = srs_dataframe.with_columns(pl.lit("SRS").alias("type"))
    combined_dataframe = pl.concat([srs_dataframe, rmd_dataframe])



    #filter combined dataframe
    combined_dataframe = (
        combined_dataframe
        .filter(
            pl.col("repeat_len") > 4
        )
        # filter out circular overlapping repeats
        .with_columns(
            (pl.col("first_repeat") + pl.col("repeat_len")).alias("left_end"),
            (pl.col("second_repeat") + pl.col("repeat_len")).alias("right_end")
        )
        .with_columns(
            pl.when(circular)
            .then(
                pl.when((pl.col("right_end")) >= seq_len)
                .then(pl.col("right_end") - seq_len)
                .otherwise(pl.col("right_end"))
            )
            .otherwise(pl.col("right_end"))
        )
        # noncircular overlapping repeats are covered by pairwise approach already
        .filter(
            ~(
                (pl.col("right_end") > pl.col("first_repeat")) &
                (pl.col("right_end") < pl.col("left_end"))
            )
        )
        .sort(["repeat_len", "first_repeat"], descending=[True, False])
        .group_by(pl.col("repeat"))
        .agg(
            pl.first("repeat_len"),
            pl.col("first_repeat"),
            pl.col("second_repeat"),
            pl.col("distance"),
            pl.first("type")
            )
        .with_columns(
            pl.concat_list(["first_repeat", "second_repeat"]).alias("positions")
        )
        # removes smaller RMDs that have the exact same positions as larger RMDs
        .sort(["repeat_len"], descending=True)
        .with_columns(
            pl.col("positions").list.eval(pl.element().unique()).alias("positions")
        )
        .group_by(pl.col("positions"))
        .agg(pl.col("repeat"),
             pl.first("repeat_len"),
             pl.first("first_repeat"),
             pl.first("second_repeat"),
             pl.first("distance"),
             pl.first("type")
             )
        .explode(pl.col("repeat"))
        .group_by(pl.col("first_repeat"))
        .head(1)

        #removes shorter versions of the same repeat that start at different positions
        .with_columns(
            pl.col("positions").list.len().alias("count")
        )
        .explode(pl.col(["first_repeat", "second_repeat", "distance"]))
        .sort(["count", "first_repeat", "repeat_len"], descending=[True, False, True])
        .group_by(pl.col("repeat"), maintain_order=True)
        .agg(
            pl.col("first_repeat"),
            pl.col("second_repeat"),
            pl.first("repeat_len"),
            pl.col("distance"),
            pl.first("count"),
            pl.col("type"),
            (pl.col("second_repeat") - seq_len).alias("negative_start")
        )

        .with_columns(
            pl.col("first_repeat").shift(1).list.unique().alias("last_first_repeat"),
            pl.col("first_repeat").list.eval(pl.element() -1).list.unique().alias("adjusted_start"),
            pl.col("repeat_len").shift(1).alias("last_len")
        )
    # delete if (last_pos == adjusted start or last_neg == adjusted start) AND last_len - len == 1
        .filter(
            (pl.col("last_first_repeat").is_null()) |
            (pl.col("adjusted_start") != pl.col("last_first_repeat")) |
            (pl.col("repeat_len") + 1 != pl.col("last_len"))
        )
    )

    filtered_df = (
    combined_dataframe
        .select(["repeat", "repeat_len", "first_repeat", "second_repeat", "distance", "type"])
        .explode(["first_repeat", "second_repeat", "distance", "type"])
        .sort(["repeat_len"], descending=True)
        .group_by(pl.col("first_repeat"), pl.col("second_repeat"))
        .head(1)
    )


    # filter out SRS nested fully inside SSRs
    # if statement needed because cross join with an empty df creates an empty df
    if ssr_dataframe.height > 0:
        ssr_ranges = (
            ssr_dataframe
            .select(
                # start is the 1st bp of SSR
                pl.col("start"),
                # end is the last bp of SSR
                (pl.col("start")+(pl.col("count")*pl.col("repeat_len")) - 1).alias("end"),
                pl.col("repeat_len"),
                pl.col("count")
            )
            .with_columns(
                pl.when(circular)
                .then(
                    pl.when(pl.col("end") >= seq_len)
                    .then(True)
                    .otherwise(False)
                    .alias("wraparound")
                )
                .otherwise(False).alias("wraparound")
            )
            .with_columns(
                (pl.col("repeat_len")*pl.col("count")).alias("ssr_length"),
                pl.when(pl.col("wraparound") == True)
                .then(pl.col("end")-seq_len)
                .otherwise(pl.col("end"))
            )
            .select("start", "end", "ssr_length", "wraparound")
        )

        # srs_dataframe with a separate row for each repeat and each SSR range
        ssr_srs_cross = filtered_df.join(ssr_ranges, how="cross")

        filtered_df = (
            ssr_srs_cross
            .with_columns(
                pl.when(pl.col("wraparound") == True)
                .then(
                    (
                        (
                            # remove if no part of srs is after end and before start of wraparound ssr
                            ~(
                                (
                                    (pl.col("first_repeat")<pl.col("start")) &
                                    ((pl.col("first_repeat")+pl.col("repeat_len")-2) > pl.col("end"))
                                ) |
                                (
                                    (pl.col("second_repeat") < pl.col("start")) &
                                    ((pl.col("second_repeat")+pl.col("repeat_len")-2) > pl.col("end"))
                                )
                            )
                        )
                        .alias("nested")
                    )
                )
                .otherwise(
                    ((pl.col("first_repeat") >= pl.col("start")) & ((pl.col("second_repeat") + pl.col("repeat_len") - 1) <= pl.col("end")))
                    .alias("nested")
                    )
            )
            .group_by(
                pl.col("first_repeat"),
                pl.col("second_repeat"),
                pl.col("repeat")
            )
            .agg(
                pl.col("nested"),
                pl.first("distance"),
                pl.first("type"),
                pl.first("repeat_len")
            )
            .filter(
                ~pl.col("nested").list.any()
            )
        )

    filtered_df = filtered_df.select("repeat", "repeat_len", "first_repeat", "second_repeat", "distance", "type")

    # split back into rmd_dataframe and srs_dataframe
    rmd_dataframe = filtered_df.filter(pl.col("type") == "RMD").drop("type")
    srs_dataframe = filtered_df.filter(pl.col("type") == "SRS").drop("type")

    return rmd_dataframe, srs_dataframe


# No longer necessary filters (already covered by pairwise approach)
# - Delete rows with overlapping repeats and only 2 occurrences (not a real repeat)
