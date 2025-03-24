import streamlit as st
import pathlib
import os
from .ingest.parse_inputs import validate_sequences
from .constants import MAX_SIZE
import polars as pl
from .webapp.SequenceState import SequenceState

class StateMachine:
    """Class for recording user session state between streamlit interactions to prevent rerunning analysis and make
    selected predictions persist."""
    def __init__(self):
        self.user_sequences = {}
        self.named_sequences = {}
        self.sequencestates = {}

    def import_sequences(self, sequences, webapp = False):
        """Import newly uploaded sequences while retaining state of existing sequences"""
        # Import sequences without overwriting old ones
        new = {seq._originhash: seq for seq in sequences}
        for key in new:
            if key in self.user_sequences:
                new[key] = self.user_sequences[key]
        if new == self.user_sequences:
            return
        self.user_sequences = new

        # Validate sequences if they changed
        validate_sequences(self.user_sequences.values(), max_len=MAX_SIZE)

        # Make webapp states
        if webapp:
            self.sequencestates = {key: SequenceState(value) for key, value in self.user_sequences.items()}

        # Update sequence names
        self.named_sequences = {}
        for i, seqhash in enumerate(self.user_sequences):
            seq = self.user_sequences[seqhash]
            if seq.description:
                sequence_name = f"{i+1}_{seq.description}"
            else:
                sequence_name = f"{i+1}_Sequence"
            self.named_sequences[sequence_name] = seqhash



    def save_results(self, folderpath, prediction_style = None, filetype = "parquet"):
        for seqname in self.named_sequences:
            seqhash = self.named_sequences[seqname]
            seqobj = self.user_sequences[seqhash]
            if not seqobj.predicted:
                if not prediction_style:
                    raise ValueError("Must specify prediction style to save results")
                elif prediction_style not in ['linear', 'pairwise']:
                    raise ValueError("Invalid prediction style: linear or pairwise")
                seqobj.call_predictions(prediction_style)
            top = seqobj.top.select(pl.exclude("predid"))
            ssrs = seqobj.ssrs.select(pl.exclude(["predid", "annotationobjects"]))
            srss = seqobj.srss.select(pl.exclude(["predid", "annotationobjects"]))
            rmds = seqobj.rmds.select(pl.exclude(["predid", "annotationobjects"]))

            folder = os.path.join(folderpath, f"{seqname}")
            path = pathlib.Path(folder)
            path.mkdir(parents=True)
            if filetype == "parquet":
                top.write_parquet(os.path.join(folder, "top.parquet"))
                ssrs.write_parquet(os.path.join(folder, "ssrs.parquet"))
                srss.write_parquet(os.path.join(folder, "srss.parquet"))
                rmds.write_parquet(os.path.join(folder, "rmds.parquet"))
            elif filetype == "csv":
                top.select(pl.exclude("annotations")).write_csv(os.path.join(folder, "top.csv"))
                ssrs.select(pl.exclude("annotations")).write_csv(os.path.join(folder, "ssrs.csv"))
                srss.select(pl.exclude("annotations")).write_csv(os.path.join(folder, "srss.csv"))
                rmds.select(pl.exclude("annotations")).write_csv(os.path.join(folder, "rmds.csv"))
            else:
                raise ValueError("Invalid filetype")
