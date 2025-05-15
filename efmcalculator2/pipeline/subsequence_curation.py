from ..constants import MIN_SHORT_SEQ_LEN, MAX_SHORT_SEQ_LEN
import polars as pl

def collect_subsequences(seq, isCircular, window_max=16) -> pl.LazyFrame:
    """Scans across a given input sequence and returns a list of subsequences

    input:
        seq (string): the sequence to be scanned
        isCircular (boolean): Whether the sequence is circular or not

    returns:
        polars dataframe containing every repeat smaller than 16 base pairs in the input sequence

    """

    seq_len = len(seq)

    if isCircular:
        # adds first 20 bp to end
        seq = seq + seq[0:20]

    def scan_genome():
        # Probably room for optimizations here
        for i, _ in enumerate(seq[:seq_len]):
            for j in range(MIN_SHORT_SEQ_LEN, MAX_SHORT_SEQ_LEN):
                if len(seq[i : i + j]) > MIN_SHORT_SEQ_LEN:
                    sub_seq = seq[i : i + j]
                    yield {"repeat": str(sub_seq), "position": i}

    repeats = (
        pl.LazyFrame(scan_genome()).group_by("repeat").agg(pl.col("position")).collect()
    ).cast({"position": pl.List(pl.Int32)})

    return repeats


def highly_mut(df: pl.DataFrame):
    df = (
        df.filter(pl.col("repeat_len") == MAX_SHORT_SEQ_LEN-1)
        .group_by(pl.col("repeat"))
        .agg(
            pl.col("pairings"), 
            pl.first("repeat_len")
            )
        .with_columns(pl.col("pairings").list.unique().alias("positions"))
        .filter(pl.col("positions").list.len() >= 40)
    )

    # terminate if there are 40+ occurrences of a single repeat
    if df.is_empty():
        return False
    else:
        return True




def _scan_RMD(df: pl.DataFrame, seq, seq_len, isCircular) -> pl.DataFrame:
    """Scans for RMDs

    input:
        df (dataframe): dataframe containing all repeats smaller than 16 base pairs
        seq (string): the sequence to be scanned
        seq_len (int): the length of the sequence
        isCircular (boolean): whether the sequence is circular or not

    returns:
        polars dataframe containing every repeat in the input sequence

    """

    known_long_repeats = df.filter(pl.col("repeat_len") == (MAX_SHORT_SEQ_LEN - 1))
    RMD_df = pl.DataFrame(
        {
            "repeat": pl.Series("repeat", [], pl.Utf8),
            "pairings": pl.Series("pairings", [], pl.List(pl.Int64)),
            "repeat_len": pl.Series("repeat_len", [], pl.Int64),
        }
    )

    def check_larger_repeats(positions, seq):
        completed = False
        step = 50
        length = MAX_SHORT_SEQ_LEN
        pos1 = positions[0]
        pos2 = positions[1]
        largest = False
        wrap = 0

        while not completed:
            prvlength = length
            length += step
            # already know they are 15 bp repeats
            if prvlength < 16:
                prvlength = 16

            # pos2 is always after pos1
            if pos2 + length > seq_len:
                if isCircular:
                    wrap += step
                    seq = seq + seq[wrap - step : wrap]
                else:
                    length = seq_len - pos2
                    largest = True
            # largest possible repeat is seq_len/2 bp
            if length >= int(seq_len / 2):
                largest = True
                length = int(seq_len / 2)
            # prevent out of bound error
            if len(seq) - pos2 < length:
                length = len(seq) - pos2

            # if sequences not equal, then largest repeat has been passed
            if (seq[pos1 : pos1 + length] != seq[pos2 : pos2 + length]) or (largest):
                # iterate 1 by 1.
                # Uses length+1 because in range stops before last index
                for j in range(prvlength, length + 1):
                    # uses j-1 because substrings end before last index
                    if seq[pos1 + (j - 1)] == seq[pos2 + (j - 1)]:
                        sub_seq = seq[pos1 : pos1 + j]
                        yield {
                            "repeat": str(sub_seq),
                            "pairings": [pos1, pos2],
                            "repeat_len": j,
                        }
                    else:
                        break
                completed = True

    def store_RMD(positions, seq):
        nonlocal RMD_df
        repeats = pl.DataFrame(check_larger_repeats(positions, seq))

        # if larger repeats were found, then repeats will have 3 columns ("repeat", "pairings", "repeat_len")
        if repeats.width == 3:
            RMD_df = pl.concat([RMD_df, repeats])
        return None

    # Apply the function to the DataFrame
    known_long_repeats.with_columns(
        pl.col("pairings").map_elements(
            lambda pairings: store_RMD(pairings, seq), return_dtype=pl.List(pl.Null)
        )
    )
    RMD_df = RMD_df.with_columns(
        pl.col("pairings").cast(pl.List(pl.Int32)),
        pl.col("repeat_len").cast(pl.Int32))

    df = pl.concat([df, RMD_df])

    return df  # In the same format as df alongside the original data
