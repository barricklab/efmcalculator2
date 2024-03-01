import argparse
import time
import pkg_resources
import logging
import Bio.SeqIO as SeqIO
import pandas as pd
import pathlib

from .short_seq_finder import predict_RMDs

from Bio.SeqRecord import SeqRecord
from typing import Union, List

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())


def _main():
    """Runs EFM calculator using the command line"""
    start_time = time.time()

    # Parse args

    parser = argparse.ArgumentParser(description="Find Short Sequences from Fasta File")
    parser.add_argument(
        "-i", "--input",
        action="store",
        dest="inpath",
        required=True,
        help="the path to fasta file",
    )
    parser.add_argument(
        "-o", "--output",
        action="store",
        dest="outpath",
        required=True,
        help="path to out csv file name",
    )
    parser.add_argument(
        "-c", "--circular",
        dest="circular",
        action="store_true",
        required=False,
        help="Is circular?",
    )
    parser.add_argument(
        "--debug",
        action = "store_true",
        dest="verbose",
        required=False,
        help="--debug/--no-debug prints debug information",
    )

    try:
        version = pkg_resources.get_distribution("efmcalculator").version
    except pkg_resources.DistributionNotFound:
        version = "unknown"

    # Sanity check filepaths

    args = parser.parse_args()

    if not pathlib.Path(args.inpath).exists():
        raise ValueError("File {} does not exist.".format(args.inpath))

    # Check to see if output is a directory
    if pathlib.Path(args.outpath).is_dir():
        raise ValueError("File {} is a directory.".format(args.outpath))

    # Set up logger
    logging.basicConfig()
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)


    logger.info("EFM Calculator version: {}".format(version))

    # Set up circular

    args.isCirc = args.circular

    # Run
    
    sequences = SeqIO.parse(args.inpath, "fasta")

    # Unpack sequences into list
    sequences = list(sequences)

    df = efmcalculator(
        sequences=sequences,
        isCircular=args.isCirc)

    df.sort_values(["Size"], ascending=[True]).to_csv(args.outpath, index=False)

    # Logging
    logger.debug(df.head())
    t_sec = round(time.time() - start_time)
    t_min, t_sec = divmod(t_sec, 60)
    t_hour, t_min = divmod(t_min, 60)
    logger.debug("\nTime taken: {}hour:{}min:{}sec".format(t_hour, t_min, t_sec))


def efmcalculator(sequences: Union[SeqRecord, List[SeqRecord]], isCircular: bool) -> pd.DataFrame:
    """Runs EFM calculator on input SeqRecords. Returns a pandas DataFrame of repeats.
    
    Args:
        sequences (Union[SeqRecord, List[SeqRecord]]): Input SeqRecords
        isCircular (bool): Treat input sequences as circular

    Returns:
        pd.DataFrame: Pandas DataFrame of repeats
    """

    data = {"Sequence": [], "Size": [], "Distance": [], "Start-Pos": [], "End-Pos": []}
    df = pd.DataFrame(data)

    if isinstance(sequences, SeqRecord):
        sequences = [sequences]

    logger.debug("Number of sequences: {}".format(len(sequences)))

    for record in sequences:
        # seq_len = len(record.seq)
        # adds first 20 bp to end
        if isCircular:
            record = record + record[0:20]
        # Strips of new line special character
        seq = record.seq.strip("\n")
        seq_len = len(seq)
        predict_RMDs(seq, df, seq_len, isCircular)
    return df

if __name__ == "__main__":
    _main()

# EOF
