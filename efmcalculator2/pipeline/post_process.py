from .mutation_rates import ssr_mut_rate_vector, srs_mut_rate_vector, rmd_mut_rate_vector
from .filtering import filter_ssrs, filter_direct_repeats, rename_tandem_repeat
from .features import assign_features_ssr, assign_features_rmd
from ..constants import THRESHOLD

import polars as pl

def post_process(ssr_df, srs_df, rmd_df, seqobj, isCircular):
    if srs_df.height == 1:
        if srs_df[0]["repeat"][0][0] == "too mutagenic":
            ssr_df = ssr_df.with_columns([pl.col("repeat").list.join(", ").alias("repeat"), pl.lit(-1.0).alias("mutation_rate")])
            srs_df = srs_df.with_columns([pl.col("repeat").list.join(", ").alias("repeat"), pl.lit(-1.0).alias("mutation_rate")])
            rmd_df = rmd_df.with_columns([pl.col("repeat").list.join(", ").alias("repeat"), pl.lit(-1.0).alias("mutation_rate")])
        else:
        # Perform Filtering
            ssr_df = filter_ssrs(ssr_df, len(seqobj), isCircular)
            rmd_df, srs_df = filter_direct_repeats(rmd_df, srs_df, len(seqobj), ssr_df, isCircular)

        # Calculate Mutation Rates
            ssr_df = ssr_mut_rate_vector(ssr_df)
            srs_df = srs_mut_rate_vector(srs_df)
            rmd_df = rmd_mut_rate_vector(rmd_df)
    else:
    # Perform Filtering
        ssr_df = filter_ssrs(ssr_df, len(seqobj), isCircular)
        rmd_df, srs_df = filter_direct_repeats(rmd_df, srs_df, len(seqobj), ssr_df, isCircular)

    # Calculate Mutation Rates
        ssr_df = ssr_mut_rate_vector(ssr_df)
        srs_df = srs_mut_rate_vector(srs_df)
        rmd_df = rmd_mut_rate_vector(rmd_df)
        

    # Filter on minimum threshold
    #ssr_df = ssr_df.filter(pl.col("mutation_rate") > THRESHOLD)
    #srs_df = srs_df.filter(pl.col("mutation_rate") > THRESHOLD)
    #rmd_df = rmd_df.filter(pl.col("mutation_rate") > THRESHOLD)


    # Classify SRS and RMD as "Tandem Repeat" when appropriate
    srs_df = rename_tandem_repeat(srs_df, seqobj, isCircular)
    rmd_df = rename_tandem_repeat(rmd_df, seqobj, isCircular)

    # Apply annotations
    if seqobj.annotations:
        ssr_df = assign_features_ssr(ssr_df, seqobj, isCircular)
        srs_df = assign_features_rmd(srs_df, seqobj, isCircular)
        rmd_df = assign_features_rmd(rmd_df, seqobj, isCircular)

    return ssr_df, srs_df, rmd_df
