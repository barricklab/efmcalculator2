import streamlit as st

class StateMachine:
    """Class for recording user session state between streamlit interactions to prevent rerunning analysis and make
    selected predictions persist."""
    def __init__(self):
        self.user_sequences = []
        self.seqobjects = []

    def check_differences(self, sequences):
        """Checks to determine whether any files have changed"""
        if not self.user_sequences:
            return True
        for seqobject in sequences:
            if seqobject not in self.seqobjects:
                return True
        for seqobject in self.seqobjects:
            if seqobject not in sequences:
                return True
        return False

    def import_sequences(self, sequences):
        checked_differences = self.check_differences(sequences)
        if not checked_differences:
            return

        # Forget sequences that have been removed
        retained_sequences = zip(self.user_sequences, self.seqobjects)
        retained_sequences = [seq for seq in retained_sequences if seq[1] in sequences]
        self.user_sequences = [seq[0] for seq in retained_sequences]
        self.seqobjects = [seq[1] for seq in retained_sequences]

        with st.spinner("Calculating..."):
            # Calculate for every new sequence
            for seq in [seq for seq in sequences if seq not in self.user_sequences]:
                pass

    def download_results(self):
        pass


class UserSequence:
    """Assistant class of state machine that handles sequence specific state data"""
    def __init__(self, seqobj):
        self.seqobj = seqobj
        self.seqhash = hash(str(seqobj.seq.strip("\n\n").upper().replace("U", "T")))
        self.name = f"{self.seqhash}_Sequence"

        self.ssr_df = None
        self.srs_df = None
        self.rmd_df = None

    def run_predictions(self):
        pass

    def assign_top_hits(self):
        pass
