import polars as pl

def filter_ssrs(ssr_dataframe, seq_len, circular):
    if len(ssr_dataframe) > 0:
        ssr_dataframe = (
            ssr_dataframe.lazy()
            # Filter based on SSR definition
            .filter(
                (pl.col("repeat_len") >= 3)
                .and_(pl.col("count") >= 3)
                .or_((pl.col("repeat_len") <= 2).and_(pl.col("count") >= 4))
            )   #! EFM1 is bugged with repeat_ len>= 3 and <= 2 bp
        ).collect()

        ssr_dataframe = (
        ssr_dataframe.sort([pl.col("start"), pl.col("repeat_len")], descending=[False, True])
            .with_columns(
                # end is the last base pair of the SSR
                (pl.col("start") + (pl.col("repeat_len")*pl.col("count"))-1).alias("end")
            )
            .with_columns(
                # if circular, creates modified start and stop positions
                pl.when(circular)
                .then(
                    pl.when((pl.col("end")) >= seq_len)
                    .then((pl.col("start")-seq_len).alias("modified_start"))
                    .otherwise(pl.col("start").alias("modified_start"))
                    )
                .otherwise(pl.col("start").alias("modified_start")),
                pl.when(circular)
                .then(
                    pl.when((pl.col("end")) >= seq_len)
                    .then((pl.col("end")-seq_len).alias("modified_end"))
                    .otherwise(pl.col("end").alias("modified_end"))
                    )
                .otherwise(pl.col("end").alias("modified_end"))
            )
            .with_row_index()
            )
        
        joined = ssr_dataframe.join(ssr_dataframe, how="cross")

        # remove rows compared to itself
        joined = joined.filter(pl.col("index") != pl.col("index_right"))
        if(not circular):
            joined = (
            joined.filter(
                    # these repeats are fully contained in another repeat, and has less count
                    (
                        (pl.col("start") >= pl.col("start_right")) &
                        (pl.col("end") <= pl.col("end_right")) &
                        (pl.col("count") <= pl.col("count_right"))
                    ) |
                    # these repeats are alternate versions of other SSRs
                    (
                        (pl.col("start") - pl.col("start_right") == pl.col("end") - pl.col("end_right")) &
                        ((pl.col("start") - pl.col("start_right")).abs() < pl.col("repeat_len")) &
                        # want to keep the first version 
                        (pl.col("start") > pl.col("start_right"))
                    )
                )
            )
        else:
            joined = (
            joined.filter(
                    # these repeats are fully contained in another repeat, and has less count
                    (
                        (
                            # use either modified start or modified end for repeats that wrap around
                            ((pl.col("start") >= pl.col("start_right")) &
                            (pl.col("end") <= pl.col("end_right"))
                            ) |
                            ((pl.col("start") >= pl.col("modified_start_right")) &
                            (pl.col("end") <= pl.col("modified_end_right"))
                            ) |
                            ((pl.col("modified_start") >= pl.col("start_right")) &
                            (pl.col("modified_end") <= pl.col("end_right"))
                            ) |
                            ((pl.col("modified_start") >= pl.col("modified_start_right")) &
                            (pl.col("modified_end") <= pl.col("modified_end_right"))
                            )
                        ) &
                        (pl.col("count") <= pl.col("count_right"))
                    ) |
                    # these repeats are alternate versions of other SSRs
                    (
                        (pl.col("start") - pl.col("start_right") == pl.col("end") - pl.col("end_right")) &
                        ((pl.col("start") - pl.col("start_right")).abs() < pl.col("repeat_len")) &
                        # want to keep the first version 
                        (pl.col("start") > pl.col("start_right"))
                    )
                )
            )

        removed_indices = joined["index"].unique()

        ssr_dataframe = (
            ssr_dataframe
            .filter(
                ~(pl.col("index").is_in(removed_indices))
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
        .sort(["first_repeat", "repeat_len"], descending=[False, True])
        .group_by(pl.col("first_repeat"), pl.col("repeat"))
        .agg(
            pl.col("second_repeat"),
            pl.first("repeat_len"),
            pl.col("distance"),
            pl.col("type")
        )
        .sort(["first_repeat", "repeat_len"], descending=[False, True])
        .with_columns(
            pl.col("first_repeat").shift(1).alias("last_first_repeat"),
            pl.col("second_repeat").shift(1).list.unique().alias("last_second_repeat"),
            pl.col("repeat_len").shift(1).alias("last_len")
        )
        .explode(["second_repeat", "distance", "type"])
        .with_row_index()
    )


    # Create a list of indices that should be deleted
    # Delete instances of smaller repeats that are in bigger repeats (ex. 2 copies of "AAGTCAT" and 3 copies of "AAGTCA". Delete the instance of "AAGTCA" that 
    # corresponds to the "AAGTCAT" repeat and keep the other 2 pairwise repeats)
    filter_out = (
        combined_dataframe
        .filter(
            (pl.col("first_repeat") == pl.col("last_first_repeat")) &
            (pl.col("last_second_repeat").list.contains(pl.col("second_repeat"))) &
            (pl.col("repeat_len") < pl.col("last_len"))
            )
        .select("index").to_series().to_list()
    )
    

    # Create another list of indices that should be deleted
    # Delete shorter versions of the same repeat that start at different positions
    filter_out_2 = (
        combined_dataframe
        .group_by("first_repeat", maintain_order=True)
        .agg(
            pl.col("repeat"),
            pl.col("second_repeat"),
            pl.first("repeat_len"),
            pl.col("distance"),
            pl.col("type"), 
            pl.first("last_first_repeat"), 
            pl.first("last_second_repeat"),
            pl.first("last_len"), 
            pl.col("index")
        )
        .with_columns(
            pl.col("first_repeat").shift(1).alias("last_first_repeat"),
            pl.col("second_repeat").shift(1).list.unique().alias("last_second_repeat"),
            pl.col("repeat_len").shift(1).alias("last_len")
        )
        .explode(["repeat", "second_repeat", "distance", "type", "index"])
        .with_columns(
            (pl.col("first_repeat") - pl.col("last_first_repeat")).alias("difference")
        )
        .with_columns(
            (pl.col("first_repeat") - pl.col("difference")).alias("adjusted_first_repeat"),
            (pl.col("second_repeat") - pl.col("difference")).alias("adjusted_second_repeat")
        )
        .filter(
            (pl.col("adjusted_first_repeat") == pl.col("last_first_repeat")) &
            (pl.col("last_second_repeat").list.contains(pl.col("adjusted_second_repeat"))) &
            (pl.col("repeat_len") <= pl.col("last_len"))
        )
        .select("index").to_series().to_list()
        )

    # Filter out repetas with indices in either of the lists
    filtered_df = combined_dataframe.filter(
        ~(
            (pl.col("index").is_in(filter_out)) |
            (pl.col("index").is_in(filter_out_2))
        )
    )


    # filter out SRS nested fully inside SSRs
    # if statement needed because cross join with an empty df creates an empty df
    if ssr_dataframe.height > 0:
        ssr_dataframe = (
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

        ssr_ranges = list(ssr_dataframe.select([pl.col("start"), pl.col("end"), pl.col("wraparound")]).iter_rows())
        

        # Function to check 'nested' condition for a single row against all tuples
        def check_nested(row):
            first_repeat = row["first_repeat"]
            second_repeat = row["second_repeat"]
            repeat_len = row["repeat_len"]

            nested_values = []
            for start, end, wraparound in ssr_ranges:
                if wraparound:
                    condition = not (
                        ((first_repeat < start) & ((first_repeat + repeat_len - 2) > end)) |
                        ((second_repeat < start) & ((second_repeat + repeat_len - 2) > end))
                    )
                else:
                    condition = (first_repeat >= start) & ((second_repeat + repeat_len - 1) <= end)
                
                nested_values.append(condition)

            return any(nested_values)  # If any tuple results in True, set nested=True


        # Apply function to each row
        filtered_df = (
            filtered_df.with_columns(
                pl.struct(["first_repeat", "second_repeat", "repeat_len"])
                .map_elements(check_nested, return_dtype=pl.Boolean)
                .alias("nested")
            )
            .filter(~pl.col("nested"))  # Filter out rows where nested=True
        )

    filtered_df = filtered_df.select("repeat", "repeat_len", "first_repeat", "second_repeat", "distance", "type")

    # split back into rmd_dataframe and srs_dataframe
    rmd_dataframe = filtered_df.filter(pl.col("type") == "RMD").drop("type")
    srs_dataframe = filtered_df.filter(pl.col("type") == "SRS").drop("type")

    return rmd_dataframe, srs_dataframe


# No longer necessary filters (already covered by pairwise approach)
# - Delete rows with overlapping repeats and only 2 occurrences (not a real repeat)
