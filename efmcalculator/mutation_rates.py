import polars as pl
from .constants import SUB_RATE


def ssr_mut_rate(repeat_count, unit_length, org):
    """
    Calculates mutation rate for simple sequence repeats
    :param repeat_count: Number of times the repeating unit occurs
    :param unit_length: Length of repeating unit
    :param org: Host organism
    :return: Mutation rate
    """
    mut_rate = float(0)
    if org == "ecoli" or org == "reca":
        if unit_length == 1:
            # Formula based on analysis of Lee et. al. data
            mut_rate = float(10 ** (0.72896 * repeat_count - 12.91471))
        elif unit_length > 1:
            mut_rate = float(10 ** (0.06282 * repeat_count - 4.74882))
    elif org == "yeast":
        if unit_length == 1:
            mut_rate = float(10 ** (0.3092 * repeat_count - 7.3220))
        elif unit_length > 1:
            mut_rate = float(10 ** (0.11141 * repeat_count - 7.65810))
    return mut_rate


def ssr_mut_rate_vector(ssr_df, org="ecoli"):
    if org == "ecoli" or org == "reca":
        ssr_df = (
            ssr_df.lazy().with_columns(
                mutation_rate=(
                    pl.when(pl.col("repeat_len") == 1)
                    .then(10 ** (0.72896 * pl.col("count") - 12.91471))
                    .otherwise(10 ** (0.06282 * pl.col("count") - 4.74882))
                )
            )
        ).collect()

    elif org == "yeast":
        ssr_df = (
            ssr_df.lazy().with_columns(
                mutation_rate=(
                    pl.when(pl.col("repeat_len") == 1)
                    .then(10 ** (0.3092 * pl.col("count") - 7.3220))
                    .otherwise(10 ** (0.11141 * pl.col("count") - 7.65810))
                )
            )
        ).collect()
    return ssr_df


def rmd_mut_rate(length, distance, org):
    """
    Calculate the recombination rate based on the Oliviera, et. al. formula
    :param length: Length of homologous region
    :param distance: Bases between the end of the first repeat and beginning of second
    :param org: Host organism
    :return: Recombination rate
    """
    spacer = distance
    # If the homologous sequences overlap we can't calculate a rate
    if spacer < 0:
        return 0
    if org == "ecoli" or org == "yeast":
        recombo_rate = float(
            ((8.8 + spacer) ** (-29.0 / length)) * (length / (1 + 1465.6 * length))
        )
    elif org == "reca":
        recombo_rate = float(
            ((200.4 + spacer) ** (-8.8 / length))
            * (length / (1 + 2163.0 * length + 14438.6 * spacer))
        )
    else:
        raise ValueError("Invalid org")

    return recombo_rate


def rmd_mut_rate_vector(rmd_df, org="ecoli"):
    if org == "ecoli" or org == "yeast":
        rmd_df = (
            rmd_df.lazy().with_columns(
                mutation_rate=(
                    ((8.8 + pl.col("distance")) ** (-29.0 / pl.col("repeat_len")))
                    * (pl.col("repeat_len") / (1 + 1465.6 * pl.col("repeat_len")))
                )
            )
        ).collect()

    elif org == "reca":
        rmd_df = (
            rmd_df.lazy().with_columns(
                mutation_rate=(
                    (200.4 + pl.col("distance")) ** (-8.8 / pl.col("repeat_len"))
                )
                * (
                    pl.col("repeat_len")
                    / (1 + 2163.0 * pl.col("repeat_len") + 14438.6 * pl.col("distance"))
                )
            )
        ).collect()
    return rmd_df


def rip_score(ssr_df, srs_df, rmd_df, sequence_length):
    if isinstance(ssr_df, pl.DataFrame) and not ssr_df.is_empty():
        ssr_sum = ssr_df.select(pl.sum("mutation_rate")).item()
    else:
        ssr_sum = 0
    if isinstance(srs_df, pl.DataFrame) and not srs_df.is_empty():
        srs_sum = srs_df.select(pl.sum("mutation_rate")).item()
    else:
        srs_sum = 0
    if isinstance(rmd_df, pl.DataFrame) and not rmd_df.is_empty():
        rmd_sum = rmd_df.select(pl.sum("mutation_rate")).item()
    else:
        rmd_sum = 0

    base_rate = float(sequence_length) * float(SUB_RATE)
    # Add in the mutation rate of an individual nucleotide
    r_sum = float(ssr_sum + srs_sum + rmd_sum + base_rate)

    # Set the maximum rate sum to 1 for now.
    if r_sum > 1:
        r_sum = float(1)
    rel_rate = float(r_sum) / float(base_rate)

    return {
        "rip": rel_rate,
        "ssr_sum": ssr_sum,
        "srs_sum": srs_sum,
        "rmd_sum": rmd_sum,
        "bps_sum": base_rate,
    }
