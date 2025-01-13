import streamlit as st
from ..parse_inputs import validate_sequences
from ..constants import MAX_SIZE

class StateMachine:
    """Class for recording user session state between streamlit interactions to prevent rerunning analysis and make
    selected predictions persist."""
    def __init__(self):
        self.user_sequences = {}
        self.named_sequences = {}

    def import_sequences(self, sequences):
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

        # Update sequence names
        self.named_sequences = {}
        for i, seqhash in enumerate(self.user_sequences):
            seq = self.user_sequences[seqhash]
            if seq.description:
                sequence_name = f"{i+1}_{seq.description}"
            else:
                sequence_name = f"{i+1}_Sequence"
            self.named_sequences[sequence_name] = seqhash

    def download_results(self):
        pass
