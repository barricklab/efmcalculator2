import argparse
import time
import logging
import Bio.SeqIO as SeqIO
import pandas as pd
import polars as pl
import pathlib

from .short_seq_finder import predict
from .SRS_filter import filter_redundant
from .filtering import filter_ssrs, filter_rmds
from .mutation_rates import ssr_mut_rate_vector, rmd_mut_rate_vector

from Bio.SeqRecord import SeqRecord
from typing import Union, List
from importlib.metadata import version, PackageNotFoundError
from rich.logging import RichHandler


FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)
logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())


def _main():
    """Runs EFM calculator using the command line"""
    start_time = time.time()

    # Parse args

    parser = argparse.ArgumentParser(description="Find Short Sequences from Fasta File")
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="inpath",
        required=True,
        help="the path to fasta/genbank file",
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="outpath",
        required=True,
        help="path to out csv file name",
    )
    parser.add_argument(
        "-s",
        "--strategy",
        action="store",
        dest="strategy",
        required=False,
        default="linear",
        help="the strategy to use for predicting deletions. Must be one of 'pairwise' or 'linear' (default: 'linear')",
    )
    parser.add_argument(
        "-c",
        "--circular",
        dest="circular",
        action="store_true",
        required=False,
        help="Is circular?",
    )
    parser.add_argument(
        "-j",
        "--threads",
        dest="threads",
        action="store",
        required=False,
        type=int,
        help="Number of threads to use",
    )
    parser.add_argument(
        "-v",
        "--verbosity",
        action="store",
        type=int,
        dest="verbose",
        required=False,
        default=1,
        help="0 - Silent | 1 Basic Information | 2 Debug",
    )

    try:
        pkgversion = version("efmcalculator")
    except PackageNotFoundError:
        pkgversion = "unknown"

    # Sanity check filepaths

    args = parser.parse_args()

    # Check to see if output is a directory
    if pathlib.Path(args.outpath).is_dir():
        raise ValueError("File {} is a directory.".format(args.outpath))

    # Set up logger
    logging.basicConfig()
    if args.verbose == 0:
        logging.root.setLevel(logging.ERROR)
    elif args.verbose == 1:
        logging.root.setLevel(logging.INFO)
    elif args.verbose == 2:
        logging.root.setLevel(logging.DEBUG)
    else:
        logger.error(f"Invalid value for verbosity flag '{args.verbose}'")
        parser.print_help()
        exit(1)

    if args.verbose < 1:
        logger.error(f"Invalid value for thread flag '{args.threads}'")
        parser.print_help()
        exit(1)

    if args.strategy not in ["pairwise", "linear"]:
        logger.error(f"Invalid value for strategy flag '{args.strategy}'")
        parser.print_help()
        exit(1)

    logger.info("EFM Calculator {}".format(pkgversion))

    # Set up circular

    args.isCirc = args.circular

    # Grab sequence information
    inpath = pathlib.Path(args.inpath)
    if not inpath.exists():
        raise ValueError("File {} does not exist.".format(args.inpath))
    elif inpath.suffix in [".fasta", ".fa"]:
        sequences = SeqIO.parse(args.inpath, "fasta")
    elif inpath.suffix in [".gb", ".gbk", ".gbff"]:
        sequences = SeqIO.parse(args.inpath, "genbank")
    else:
        raise ValueError(
            "File {} is not a known file format. Must be in FASTA or GenBank.".format(
                args.inpath
            )
        )

    # Unpack sequences into list
    sequences = list(sequences)

    df = efmcalculator(
        sequences=sequences, strategy=args.strategy, isCircular=args.isCirc
    )
    df = filter_redundant(df)

    # Index every position in the dataframe up by 1 @TODO

    # Output
    df.to_csv(args.outpath, index=False)

    # Logging
    logger.debug(df.head())
    t_sec = round(time.time() - start_time)
    t_min, t_sec = divmod(t_sec, 60)
    t_hour, t_min = divmod(t_min, 60)
    logger.info(f"EFMCalculator completed in {t_hour}:{t_min}:{t_sec}")


def efmcalculator(
    sequences: Union[SeqRecord, List[SeqRecord]],
    strategy: str,
    isCircular: bool,
    threads: int = 1,
) -> pd.DataFrame:
    """Runs EFM calculator on input SeqRecords. Returns a pandas DataFrame of repeats.

    Args:
        sequences (Union[SeqRecord, List[SeqRecord]]): Input SeqRecords
        strategy (str): The strategy to use for predicting deletions. Must be one of "pairwise" or "linear".
        isCircular (bool): Treat input sequences as circular
        threads: Number of threads to use

    Returns:
        pd.DataFrame: Pandas DataFrame of repeats
    """

    data = {
        "Sequence": [],
        "Size": [],
        "Distance": [],
        "Start-Pos": [],
        "End-Pos": [],
        "Classifier": [],
        "Repeat Rate": [],
        "Mutation Rate": [],
    }
    df = pd.DataFrame(data)

    if isinstance(sequences, SeqRecord):
        sequences = [sequences]

    for record in sequences:
        logger.info("Running on {}".format(record.name))

        # Perform predictions
        seq = record.seq.strip("\n").upper()
        ssr_df, rmd_df = predict(seq, df, strategy, isCircular, threads)

        # Perform Filtering
        ssr_df = filter_ssrs(ssr_df)
        rmd_df = filter_rmds(rmd_df)

        # Calculate Mutation Rates

        ssr_df = ssr_mut_rate_vector(ssr_df)
        rmd_df = rmd_mut_rate_vector(rmd_df)

        print("ssr")
        print(ssr_df)
        print("rmd")
        print(rmd_df)

    # -------------------------------- Legacy code to be removed

    # Create dataframe of observed repeats, rather than of observed sequences that have duplicates

    df["position"] = list(
        zip(
            df["Start-Pos"],
            df["End-Pos"],
            df["Classifier"],
            df["Repeat Rate"],
            df["Mutation Rate"],
        )
    )
    results = df.groupby("Sequence").agg("position").agg(["unique"])

    consolidated_df = pd.DataFrame(results)
    consolidated_df = consolidated_df.reset_index().rename(
        columns={"Sequence": "sequence", "unique": "positions"}
    )
    consolidated_df["length"] = consolidated_df["sequence"].apply(lambda x: len(x))
    consolidated_df["occurrences"] = consolidated_df["positions"].apply(
        lambda x: len(x)
    )

    # Drop repeats with one occurance
    RMD_df = consolidated_df["occurrences"] > 1
    SSR_df = consolidated_df["positions"].apply(
        lambda pos: any("SSR" in p for p in pos)
    )
    consolidated_df = consolidated_df[RMD_df | SSR_df]

    return consolidated_df

    # --------------------------------- End legacy code


if __name__ == "__main__":
    _main()

# EOF
