import logging
import polars as pl
import tempfile
import itertools
from progress.spinner import Spinner
from progress.bar import IncrementalBar
import numpy as np

class SeqAttr:
    def __init__(self, sub_seq, distance, start_pos, end_pos, note):
        self.sub_seq = str(sub_seq)
        self.length = len(sub_seq)
        self.distance = distance
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.note = note

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())

def pairwise_list_column(polars_df, column) -> pl.DataFrame:
    """Recieve a polars dataframe with column of [List[type]]
    Returns the dataframe back with a pairwise list of positions"""

    # Shamelessly adapted from https://stackoverflow.com/questions/77354114/how-would-i-generate-combinations-of-items-within-polars-using-the-native-expres

    nrows = polars_df.select(pl.len()).item()

    if logger.isEnabledFor(logging.INFO):
        bar = IncrementalBar('Generating candidate deletions', max=nrows)
    else:
        bar = FakeBar()

    def map_function(list_o_things):
        bar.next()
        return [sorted((thing_1, thing_2))
                for thing_1, thing_2 
                in itertools.combinations(list_o_things, 2)]

    pairwise = (polars_df
    .lazy()
    .with_columns(pl.col(column)
    .map_elements(map_function,
        return_dtype = pl.List(pl.List(pl.Int64))).alias(f'{column}_pairwise')
    )
    )

    pairwise = pairwise.collect()
    bar.finish()

    return pairwise

    # There's probably a way to optimize this but the
    # stackoverflow answer is wrong, doesn't work for >1 row

def _calculate_distances(polars_df, seq_len, circular) -> pl.DataFrame:


    distance_df = polars_df.with_columns(
        distance = pl.col('position_pairwise').list.diff().list.get(1) - pl.col('repeat_len')
        )

    if circular:
        distance_df = distance_df.with_columns(
            distance = pl.when(pl.col('distance') > seq_len / 2).then(
                        seq_len - pl.col('distance') + pl.col('repeat_len')
                        ).otherwise(pl.col('distance')),
            wraparound = pl.when(pl.col('distance') > seq_len / 2).then(
                        True
                        ).otherwise(False)
        )
    
    return distance_df

def _categorize_efm(polars_df) -> pl.DataFrame:

    categorized_df = polars_df.with_columns(
        is_RMD = pl.when(pl.col('distance') > 0).then(True).
        otherwise(False)
    )

    return categorized_df


def _collapse_ssr(polars_df) -> pl.DataFrame:


    collapsed_ssrs = (
        polars_df.filter(pl.col('is_RMD') == False)
        .select(['repeat', 'repeat_len', 'position_pairwise'])
        .lazy()
        .explode('position_pairwise')  # Collect the positions from all potentially participating SSRs
        .group_by(['repeat', 'repeat_len'])
        .agg('position_pairwise')
        .with_columns(
            positions = pl.col('position_pairwise').list
            .unique().list.sort(),
        )
        .with_columns(
            differences = pl.col('positions').list.diff()
        )
        .collect()
    )
    # Somehow couldnt figure out how to do this in pure polars

    # Identify start positions. If distance !=0, its not a start
    collapsed_ssrs = collapsed_ssrs.to_pandas()
    collapsed_ssrs['differences'] = (collapsed_ssrs['differences'] - collapsed_ssrs['repeat_len'] ).fillna(0)
    collapsed_ssrs = pl.from_pandas(collapsed_ssrs)

    collapsed_ssrs = (
        collapsed_ssrs.with_columns(
            truth_table = (pl.col('differences')
                    .list.eval((pl.element() != 0)
                    .or_(pl.element().is_null())))
        )
        )

    # Apply the truth table
    collapsed_ssrs = collapsed_ssrs.to_pandas()
    collapsed_ssrs['starts'] = (collapsed_ssrs['truth_table']) * (collapsed_ssrs['positions']+1)


    # Fill out 
    collapsed_ssrs = pl.from_pandas(collapsed_ssrs)
    collapsed_ssrs = collapsed_ssrs.with_columns(
        starts = pl.col('starts').explode().replace(0, None).forward_fill().implode().over('repeat')
    )

    # We had to apply an offeset with the truth table
    # To prevent bp-0 start positions from nulling out. This undoes that.
    collapsed_ssrs = collapsed_ssrs.to_pandas()
    collapsed_ssrs['starts'] = collapsed_ssrs['starts'] -1
    collapsed_ssrs = pl.from_pandas(collapsed_ssrs).select(pl.col(["repeat",
                                                                   "repeat_len",
                                                                    "starts"]))

    # Count repeats from each start position
    collapsed_ssrs = (
        collapsed_ssrs
        .explode('starts')
        .group_by('repeat', 'repeat_len', 'starts').count())


    # Filter based on SSR definition
    collapsed_ssrs = (
        collapsed_ssrs.filter(
            (pl.col('repeat_len') >= 2).and_(pl.col('count') >=3)
            .or_((pl.col('repeat_len') == 1).and_(pl.col('count') >=4))
        )
    )



    # @TODO: Correct for circular

    return polars_df

def _apply_mutation_rates(polars_df) -> pl.DataFrame:
    return polars_df



def _build_seq_attr(sub_seq, seq_len, start_positions, isCircular, count):
    distance = 1
    start_pos = 0
    rem_start = 0
    end_pos = 0
    prv_end_pos = 0
    prv_start_pos = 0


    for start_pos in start_positions:
        note = ""
        end_pos = start_pos + len(sub_seq)
        #if end_pos >= seq_len:
         #   end_pos = end_pos - seq_len
        start_pos += 1
        distance = start_pos - prv_end_pos

        # fixes distance for
        if isCircular == True:
            if distance > seq_len / 2:
                distance = (seq_len + prv_start_pos) - end_pos
            if start_pos < seq_len and end_pos > seq_len: # if repeat wraps around
                # fix end_pos
                end_pos = end_pos - seq_len
        if distance == 1:
            note = "SSR"
        # if overlapping
        if distance < 1:
            note = "skip for SSR"
            if count == 2:
                note = "skip"
        yield SeqAttr(sub_seq, distance, start_pos, end_pos, note)
        rem_start = start_pos + 1
        prv_end_pos = end_pos
        prv_start_pos = start_pos


def _find_repeat_positions(seq, sub_seq, seq_len, isCircular, count):
    '''testing function for sequence scanning'''
    start_pos = 0
    rem_start = 0
    end_pos = 0
    prv_start_pos = 0
    start_positions = []

    start_pos = seq.find(sub_seq, rem_start)
    while start_pos != -1:
        start_positions.append(start_pos)
        start_pos += 1
        rem_start = start_pos + 1

        if len(sub_seq) == 1:
            start_pos = seq.find(sub_seq, rem_start-1)
        else:
            start_pos = seq.find(sub_seq, rem_start)
    return start_positions


