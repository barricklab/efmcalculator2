import argparse
import time
import logging
import Bio.SeqIO as SeqIO
import pandas as pd
import polars as pl
import pathlib
import os
from Bio.Seq import Seq

from .short_seq_finder import predict
from .SRS_filter import filter_redundant
from .filtering import filter_ssrs, filter_direct_repeats
from .mutation_rates import ssr_mut_rate_vector, srs_mut_rate_vector, rmd_mut_rate_vector, rip_score
from .constants import VALID_STRATEGIES, FASTA_EXTS, GBK_EXTS
from .parse_inputs import parse_file, validate_sequences, BadSequenceError

from .utilities import is_pathname_valid, is_path_creatable, sanitize_filename
from .visualization.graph import make_plot
from .visualization.make_webpage import make_standalone_page, export_html

from Bio.SeqRecord import SeqRecord
from typing import Union, List, Set, Generator
from importlib.metadata import version, PackageNotFoundError
from rich.logging import RichHandler
from streamlit_extras.concurrency_limiter import concurrency_limiter
import streamlit as st

from bokeh.plotting import output_file, save

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

    parser = argparse.ArgumentParser(description="Find Mutation Hotspots from Input Sequences")
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="inpath",
        required=True,
        help="the path to fasta/genbank/csv file",
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="outpath",
        required=True,
        help="path to output prefix",
    )
    parser.add_argument(
        "-s",
        "--strategy",
        action="store",
        dest="strategy",
        required=False,
        default="pairwise",
        help="the strategy to use for predicting deletions. Must be one of 'pairwise' or 'linear' (default: 'pairwise')",
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
        "-v",
        "--verbosity",
        action="store",
        type=int,
        dest="verbose",
        required=False,
        default=1,
        help="0 - Silent | 1 Basic Information | 2 Debug",
    )
    parser.add_argument(
        "--no-vis",
        action="store_true",
        dest="no_vis",
        required=False,
        help="Skip visualization",
    )

    try:
        pkgversion = version("efmcalculator")
    except PackageNotFoundError:
        pkgversion = "unknown"

    logger.info("EFM Calculator {}".format(pkgversion))

    args = parser.parse_args()

    # Set up logger  -----------
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

    # Sanity checks  ------------
    if args.strategy not in VALID_STRATEGIES:
        logger.error(f"Invalid value for strategy flag '{args.strategy}'")
        parser.print_help()
        exit(1)

    if not is_pathname_valid(args.inpath):
        logger.error(f"File {args.inpath} is not a valid path.")
        exit(1)
    elif not is_pathname_valid(args.outpath):
        logger.error(f"File {args.outpath} is not a valid path.")
        exit(1)
    elif not is_path_creatable(args.outpath):
        logger.error(f"Cannot write to {args.outpath}")
        exit(1)

    # Set up circular ------------

    args.isCirc = args.circular

    # Grab sequence information --------
    inpath = pathlib.Path(args.inpath)
    try:
        sequences = parse_file(inpath)
    except ValueError as e:
        logger.error(e)
        exit(1)
    except OSError as e:
        logger.error(e)
        exit(1)

    try:
        validate_sequences(sequences)
    except BadSequenceError as e:
        logger.error(e)
        exit(1)

    # Unpack sequences into list ---------
    sequences = list(sequences)

    # Prepend numbers to sequence names if there are more than one ---------

    if len(sequences) > 1:
        for i, seq in enumerate(sequences):
            sequences[i].name = f"{i+1}_{seq.name}"

    # Run predictions -----------
    results_generator = predict_many(
                sequences=sequences, strategy=args.strategy, isCircular=args.isCirc
            )

    results_list = []

    for result, sequence in zip(results_generator, sequences):

        logger.info(rip_score(result[0], result[1], result[2], len(sequence.seq)))
        result.append(results_list)

    bulk_output(results_list, sequences, args.outpath, args.no_vis)

    # Logging
    t = time.time() - start_time
    t_sec = int(t)
    t_msec = round((t - t_sec) * 1000)
    t_min, t_sec = divmod(t_sec, 60)
    t_hour, t_min = divmod(t_min, 60)
    logger.info(
        f"EFMCalculator completed in {t_hour:02d}h:{t_min:02d}m:{t_sec:02d}s:{t_msec:02d}ms"
    )

def bulk_output(results_list, sequences, outpath, skip_vis = False):
    for i, result in enumerate(results_list):
        input_sequence = sequences[i]

        if len(sequences) > 1:
            # https://stackoverflow.com/questions/7406102/create-sane-safe-filename-from-any-unsafe-string
            sanitized_record_name = sanitize_filename(input_sequence.name)
            folder = outpath + "/" + sanitized_record_name + "/"
        else:
            folder = outpath + "/"

        # Create output folder if it doesn't exist ---------
        try:
            if not os.path.exists(folder):
                os.makedirs(folder)
        except FileExistsError:
            logger.error(f"Output folder {folder} already exists as a file.")
            exit(1)

        # Export results ------------

        result[0].write_csv(folder + "ssr.csv")
        result[1].write_csv(folder + "srs.csv")
        result[2].write_csv(folder + "rmd.csv")

        # Run data vis ---------
        if skip_vis:
            continue
        fig, tables = make_plot(
            input_sequence, ssr=result[0], srs=result[1], rmd=result[2]
        )
        summary = rip_score(result[0], result[1], result[2], len(input_sequence.seq))
        layout = make_standalone_page(fig, tables, summary)
        export_html(layout, f"{folder}plot.html")


@concurrency_limiter(max_concurrency=1)
def predict_many(
    sequences: List[SeqRecord],
    strategy: str,
    isCircular: bool,
) -> Generator[List[pl.DataFrame], None, None]:
    """Runs EFM calculator on input SeqRecords. Generates Dataframes (SSR, SRS, and RMD) in the same order as the input sequences.

    Args:
        sequences (List[SeqRecord]): Input SeqRecords
        strategy (str): The strategy to use for predicting deletions. Must be one of "pairwise" or "linear".
        isCircular (bool): Treat input sequences as circular

    Returns:
        Set[pd.DataFrame]: Polars dataframes of repeats
    """

    # Don't trust the user

    if len(sequences) == 0:
        raise ValueError("No sequences provided")
    elif not isinstance(strategy, str) or strategy not in VALID_STRATEGIES:
        raise ValueError(f"Invalid strategy. Must be one of {VALID_STRATEGIES}")
    elif not isinstance(sequences, List) or not all(
        isinstance(seq, SeqRecord) for seq in sequences
    ):
        raise ValueError("Input sequences must be a list of SeqRecords")
    elif not isinstance(isCircular, bool):
        raise ValueError("isCircular must be a boolean")

    # Run the EFM calculator on each sample

    for record in sequences:
        logger.info("Running on {}".format(record.name))

        # Perform predictions
        seq = str(record.seq.strip("\n").upper().replace("U", "T"))
        ssr_df, srs_df, rmd_df = predict(seq, strategy, isCircular)
        ssr_df, srs_df, rmd_df = post_process(ssr_df, srs_df, rmd_df, len(record), isCircular)

        yield [ssr_df, srs_df, rmd_df]


def post_process(ssr_df, srs_df, rmd_df, seq_len, isCircular):
    # Perform Filtering
    ssr_df = filter_ssrs(ssr_df)
    rmd_df, srs_df = filter_direct_repeats(rmd_df, srs_df, seq_len, ssr_df, isCircular)

    # Calculate Mutation Rates

    ssr_df = ssr_mut_rate_vector(ssr_df)
    srs_df = srs_mut_rate_vector(srs_df)
    rmd_df = rmd_mut_rate_vector(rmd_df)

    return ssr_df, srs_df, rmd_df


if __name__ == "__main__":
    _main()

# EOF
