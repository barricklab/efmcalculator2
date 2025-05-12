
import argparse
import logging
import pathlib
import time
from importlib.metadata import version, PackageNotFoundError
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .utilities import (
    is_path_creatable,
    is_pathname_valid,
)
from .ingest.EFMSequence import EFMSequence
from .ingest.bad_state_mitigation import BadSequenceError
from .ingest.parse_inputs import (
    parse_file,
    validate_sequence,
    validate_sequences,
)

from .constants import VALID_STRATEGIES, VALID_FILETYPES
from .StateMachine import StateMachine

import logging
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)
logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())

def main():
    """CLI entry point for EFM Calculator"""

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
        "-f",
        "--filetype",
        dest="filetype",
        action="store",
        required=False,
        help="Output filetype for tables (csv | parquet)",
    )
    parser.add_argument(
        "-j",
        "--maxthreads",
        dest="threads",
        action="store",
        type=int,
        required=False,
        help="Maximum number of threads (>0)",
    )
    parser.add_argument(
        "-t",
        "--tall",
        dest="tall",
        action="store_true",
        required=False,
        help="Parallelize across samples rather than within samples",
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
        "--summary",
        dest="summaryonly",
        action="store_true",
        required=False,
        help="Only save summary information rather than all . Useful for very tall inputs.",
    )

    try:
        from ._version import version_tuple
        pkgversion = f"{version_tuple[0]}.{version_tuple[1]}.{version_tuple[2]}"
    except:
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

    if args.filetype and args.filetype not in VALID_FILETYPES:
        logger.error(f"Invalid value for filetype flag '{args.filetype}'")
        parser.print_help()
        exit(1)
    else:
        args.filetype = "csv"

    if not is_pathname_valid(args.inpath):
        logger.error(f"File {args.inpath} is not a valid path.")
        exit(1)
    elif not is_pathname_valid(args.outpath):
        logger.error(f"File {args.outpath} is not a valid path.")
        exit(1)
    elif not is_path_creatable(args.outpath):
        logger.error(f"Cannot write to {args.outpath}")
        exit(1)


    if args.tall:
        os.environ["POLARS_MAX_THREADS"] = "1"
    if args.threads is not None and args.threads <=0:
        logger.error("Max threads must be greater than 0")
        exit(1)
    elif args.tall and not args.threads:
        threads = os.cpu_count()
    else:
        threads = args.threads
        os.environ["POLARS_MAX_THREADS"] = str(threads)

    global pl
    import polars as pl


    # Set up circular ------------

    args.isCirc = args.circular

    # Grab sequence information --------
    try:
        sequences = parse_file(pathlib.Path(args.inpath), iscircular=args.circular)
    except ValueError as e:
        logger.error(e)
        exit(1)
    except OSError as e:
        try:
            validate_sequence(EFMSequence(SeqRecord(Seq(args.inpath), id="text input", name="text input"), is_circular=args.circular))
            sequences = [EFMSequence(SeqRecord(Seq(args.inpath), id="text input", name="text input"), is_circular=args.circular)]
        except:
            logger.error("Input is not an existing file or valid sequence")
            exit(1)


    # Unpack sequences into list ---------
    sequences = list(sequences)

    # Run EFM Calculator ----------------
    statemachine = StateMachine()

    try:
        statemachine.import_sequences(sequences)
    except BadSequenceError as e:
        logger.error(e)
        exit(1)

    if args.tall:
        statemachine.predict_tall(strategy=args.strategy,
                                  outpath=args.outpath,
                                  filetype=args.filetype,
                                  threads=threads,
                                  keepmem=not args.summaryonly,
                                  summaryonly=args.summaryonly)
    else:
        for i, seqobject in enumerate(statemachine.user_sequences.values()):
            logger.info(msg=f"Running on sequence {i}: {str(seqobject)}")
            seqobject.call_predictions(strategy=args.strategy)
        statemachine.save_results(args.outpath, filetype=args.filetype, summaryonly=args.summaryonly)


    # Done ------------------------------

    # Logging
    t = time.time() - start_time
    t_sec = int(t)
    t_msec = round((t - t_sec) * 1000)
    t_min, t_sec = divmod(t_sec, 60)
    t_hour, t_min = divmod(t_min, 60)
    logger.info(
        f"EFMCalculator completed in {t_hour:02d}h:{t_min:02d}m:{t_sec:02d}s:{t_msec:02d}ms"
    )

if __name__ == "__main__":
    main()
