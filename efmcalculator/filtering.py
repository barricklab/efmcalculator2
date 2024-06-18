import polars as pl

def filter_ssrs(ssrs_dataframe):
    ssrs_dataframe = (ssrs_dataframe.lazy()

        # Filter based on SSR definition
        .filter(
            (pl.col('repeat_len') >= 2).and_(pl.col('count') >=3)
            .or_((pl.col('repeat_len') == 1).and_(pl.col('count') >=4))
        )

        # Keep SSR with lowest length @ TODO


    ).collect()
    return ssrs_dataframe


def filter_rmds(rmd_dataframe):

     # Delete redundant SRS repeats @TODO

    return rmd_dataframe


# No longer necessary filters (already covered by pairwise approach)
# - Delete rows with overlapping repeats and only 2 occurrences (not a real repeat) 
