import logging
import polars as pl
import tempfile
import itertools
from progress.spinner import Spinner
from progress.bar import IncrementalBar

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
    Returns the dataframe back with the column as all combinations"""

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
        return_dtype = pl.List(pl.List(pl.Int64)))
    )
    .collect()
    )

    bar.finish()

    return pairwise

    # There's probably a way to optimize this but the
    # stackoverflow answer is wrong, doesn't work for >1 row

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


