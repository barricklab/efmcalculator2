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
        new = {seq._originhash: seq for seq in sequences}
        for key in new:
            if key in self.user_sequences:
                new[key] = self.user_sequences[key]
        if new == self.user_sequences:
            return
        self.user_sequences = new

        validate_sequences(self.user_sequences.values(), max_len=MAX_SIZE)


    def download_results(self):
        pass
