import argparse
import time
import logging
import Bio.SeqIO as SeqIO
import pandas as pd
import polars as pl
import pathlib
from rich import print
import os
from Bio.Seq import Seq

from .short_seq_finder import predict
from .filtering import filter_ssrs, filter_direct_repeats
from .mutation_rates import ssr_mut_rate_vector, srs_mut_rate_vector, rmd_mut_rate_vector, rip_score
from .constants import VALID_STRATEGIES, FASTA_EXTS, GBK_EXTS, THRESHOLD
from .parse_inputs import parse_file, validate_sequences, validate_sequence, BadSequenceError

from .utilities import is_pathname_valid, is_path_creatable, sanitize_filename

from Bio.SeqRecord import SeqRecord
from typing import Union, List, Set, Generator
from importlib.metadata import version, PackageNotFoundError
from rich.logging import RichHandler
from streamlit_extras.concurrency_limiter import concurrency_limiter
import streamlit as st

from .features import assign_features_ssr, assign_features_rmd

from bokeh.plotting import output_file, save

FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)
logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())


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
        ssr_df, srs_df, rmd_df = post_process(ssr_df, srs_df, rmd_df, record, isCircular)

        yield [ssr_df, srs_df, rmd_df]


def post_process(ssr_df, srs_df, rmd_df, seqobj, isCircular):
    # Perform Filtering
    ssr_df = filter_ssrs(ssr_df, len(seqobj), isCircular)
    rmd_df, srs_df = filter_direct_repeats(rmd_df, srs_df, len(seqobj), ssr_df, isCircular)

    # Calculate Mutation Rates

    ssr_df = ssr_mut_rate_vector(ssr_df)
    srs_df = srs_mut_rate_vector(srs_df)
    rmd_df = rmd_mut_rate_vector(rmd_df)

    # Filter on minimum threshold
    #ssr_df = ssr_df.filter(pl.col("mutation_rate") > THRESHOLD)
    #srs_df = srs_df.filter(pl.col("mutation_rate") > THRESHOLD)
    #rmd_df = rmd_df.filter(pl.col("mutation_rate") > THRESHOLD)

    # Apply annotations
    if seqobj.annotations:
        ssr_df = assign_features_ssr(ssr_df, seqobj, isCircular)
        srs_df = assign_features_rmd(srs_df, seqobj, isCircular)
        rmd_df = assign_features_rmd(rmd_df, seqobj, isCircular)

    return ssr_df, srs_df, rmd_df



# EOF
