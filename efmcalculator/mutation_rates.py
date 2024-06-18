import polars as pl

def ssr_mut_rate(repeat_count, unit_length, org):
    '''
    Calculates mutation rate for simple sequence repeats
    :param repeat_count: Number of times the repeating unit occurs
    :param unit_length: Length of repeating unit
    :param org: Host organism
    :return: Mutation rate
    '''
    mut_rate = float(0)
    if org == 'ecoli' or org == 'reca':
        if unit_length == 1:
            # Formula based on analysis of Lee et. al. data
            mut_rate = float(10 ** (0.72896 * repeat_count - 12.91471))
        elif unit_length > 1:
            mut_rate = float(10 ** (0.06282 * repeat_count - 4.74882))
    elif org == 'yeast':
        if unit_length == 1:
            mut_rate = float(10 ** (0.3092 * repeat_count - 7.3220))
        elif unit_length > 1:
            mut_rate = float(10 ** (0.11141 * repeat_count - 7.65810))
    return mut_rate

def ssr_mut_rate_vector(ssr_df, org='ecoli'):

    if org == 'ecoli' or org == 'reca':
        ssr_df = (ssr_df.lazy().with_columns(
                    mutation_rate = (
                        pl.when(pl.col('repeat_len') == 1)
                        .then(10 ** (0.72896 * pl.col('repeat_len') - 12.91471))
                        .otherwise(10 ** (0.06282 * pl.col('repeat_len') - 4.74882))
                    )
                  )
        ).collect()

    elif org == 'yeast':
        ssr_df = (ssr_df.lazy().with_columns(
                    mutation_rate = (
                        pl.when(pl.col('repeat_len') == 1)
                        .then(10 ** (0.3092 * pl.col('repeat_len') - 7.3220))
                        .otherwise(10 ** (0.11141 * pl.col('repeat_len') - 7.65810))
                    )
                  )
        ).collect()
    return ssr_df


def rmd_mut_rate(length, distance, org):
    '''
    Calculate the recombination rate based on the Oliviera, et. al. formula
    :param length: Length of homologous region
    :param distance: Bases between the end of the first repeat and beginning of second
    :param org: Host organism
    :return: Recombination rate
    '''
    spacer = distance
    # If the homologous sequences overlap we can't calculate a rate
    if spacer < 0:
        return 0
    if org == 'ecoli' or org == 'yeast':
        recombo_rate = float(((8.8 + spacer) ** (-29.0 / length)) * (length / (1 + 1465.6 * length)))
    elif org == 'reca':
        recombo_rate = float(
            ((200.4 + spacer) ** (-8.8 / length)) * (length / (1 + 2163.0 * length + 14438.6 * spacer)))

    return recombo_rate

def rmd_mut_rate_vector(rmd_df, org='ecoli'):
    if org == 'ecoli' or org == 'yeast':
        rmd_df = (rmd_df.lazy().with_columns(
                    mutation_rate = (
                       ((8.8 + pl.col('distance')) ** (-29.0 / pl.col('repeat_len'))
                       ) * (pl.col('repeat_len') / (
                        1 + 1465.6 * pl.col('repeat_len')))
                    )
                  )
        ).collect()

    elif org == 'reca':
        rmd_df = (rmd_df.lazy().with_columns(
                    mutation_rate = (
                        (200.4 + pl.col('distance')) ** (-8.8 / pl.col('repeat_len'))) * (
                        pl.col('repeat_len') / (
                        1 + 2163.0 * pl.col('repeat_len') + 14438.6 * pl.col('distance'))
                        )
                  )
        ).collect()
    return rmd_df