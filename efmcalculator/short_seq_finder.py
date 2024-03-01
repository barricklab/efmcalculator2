from Bio import SeqIO
import logging
import pkg_resources
from Bio.Seq import Seq
import pandas as pd
import time
import numpy
import short_seq_scan as sss

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())

_start_time = time.time()
_debug = False
MIN_SHORT_SEQ_LEN = 5
MAX_SHORT_SEQ_LEN = 17


def main(filename, output, isCircular):
    _isCircular = isCircular
    data = {"Sequence": [], "Size": [], "Distance": [], "Start-Pos": [], "End-Pos": []}
    df = pd.DataFrame(data)
    seq_parser(filename, output, df, isCircular)

    # Logging
    logger.debug(df.head())
    t_sec = round(time.time() - _start_time)
    t_min, t_sec = divmod(t_sec, 60)
    t_hour, t_min = divmod(t_min, 60)
    logger.debug("\nTime taken: {}hour:{}min:{}sec".format(t_hour, t_min, t_sec))


def seq_parser(filename, output, df, isCircular):
    sequences = SeqIO.parse(filename, "fasta")
    for record in sequences:
        # seq_len = len(record.seq)
        # adds first 20 bp to end
        if isCircular:
            record = record + record[0:20]
        # Strips of new line special character
        seq = record.seq.strip("\n")
        seq_len = len(seq)
        build_sub_seq_from_seq(seq, df, seq_len, isCircular)
    df.sort_values(["Size"], ascending=[True]).to_csv(output, index=False)


def build_sub_seq_from_seq(seq, df, seq_len, isCircular):
    global visited_sequences
    visited_sequences = set()
    for i, letter in enumerate(seq):
        for j in range(MAX_SHORT_SEQ_LEN, MIN_SHORT_SEQ_LEN, -1):
            if len(seq[i : i + j]) > MIN_SHORT_SEQ_LEN:
                sub_seq = seq[i : i + j]
                rem_seq = seq[i + 1 :]
                find_short_seq(seq, sub_seq, df, seq_len, isCircular)


"""
    #circular short sequence look up
    k = MIN_SHORT_SEQ_LEN
    for j in (MIN_SHORT_SEQ_LEN, MAX_SHORT_SEQ_LEN,1):
        if (j > k):
            for l in range(len(seq) % k):
                roll_arr = numpy.roll([c for c in seq],l)
                roll_seq = ''.join(str(x) for x in roll_arr)
                sub_seq = Seq(roll_seq[0:k])
                find_short_seq(seq,sub_seq,df)
"""


def find_short_seq(seq, sub_seq, df, seq_len, isCircular):
    count = seq.count_overlap(sub_seq)
    if count > 1 and str(sub_seq) not in visited_sequences:
        df["Sequence"] = df["Sequence"].astype(str)
        # Filter for largest SRS only
        # df_filter = df[df['Sequence'].str.contains(str(sub_seq))]
        # if df_filter.empty:

        for seq_attr in sss.scan_short_sequence(
            seq, sub_seq, seq_len, isCircular, count
        ):
            logger.debug(
                "Sequence = {} : Size = {} : Distance = {} : Start-Pos = {} : End-Pos = {}".format(
                    seq_attr.sub_seq,
                    seq_attr.length,
                    seq_attr.distance,
                    seq_attr.start_pos,
                    seq_attr.end_pos,
                )
            )
            df.loc[seq_attr.sub_seq + "~" + str(seq_attr.start_pos)] = [
                seq_attr.sub_seq,
                seq_attr.length,
                seq_attr.distance,
                seq_attr.start_pos,
                seq_attr.end_pos,
            ]
            if count <= 1 or str(sub_seq) not in visited_sequences:
                visited_sequences.add(str(sub_seq))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Find Short Sequences from Fasta File")
    parser.add_argument(
        "--input_file_name",
        metavar="path",
        required=True,
        help="the path to fasta file",
    )
    parser.add_argument(
        "--output_csv_file_name",
        metavar="path",
        required=True,
        help="path to out csv file name",
    )
    parser.add_argument(
        "--isCircular",
        metavar="circular",
        default="True",
        required=False,
        help="Is circular?",
    )
    parser.add_argument(
        "--debug",
        metavar="verbose",
        required=False,
        action=argparse.BooleanOptionalAction,
        help="--debug/--no-debug prints debug information",
    )

    args = parser.parse_args()
    main(
        filename=args.input_file_name,
        output=args.output_csv_file_name,
        isCircular=args.isCircular,
    )

# Example usage:
