import streamlit as st
import subprocess
import io
import glob
import Bio.SeqIO as SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation
from Bio.Seq import Seq
import base64
import zipfile
import hmac
import polars as pl
import os
from ..short_seq_finder import predict
from ..constants import VALID_EXTS, MAX_SIZE
from ..parse_inputs import parse_file, validate_sequences
from ..bad_state_mitigation import BadSequenceError
from importlib.metadata import version
from .bokeh_plot import bokeh_plot

from ..efmcalculator import post_process
from bokeh.embed import file_html
import streamlit.components.v1 as components
from tempfile import TemporaryDirectory
from pathlib import Path
from streamlit_extras.stylable_container import stylable_container
from streamlit_extras.add_vertical_space import add_vertical_space
from ..efmcalculator import predict_many
from ..mutation_rates import rip_score
from .vis_utils import eval_top
from .state_machine import StateMachine
from ..EFMSequence import EFMSequence

import hashlib

ASSET_LOCATION = os.path.join(os.path.dirname(__file__), "../visualization/assets")

def check_password():
    """Returns `True` if the user had the correct password."""

    def password_entered():
        """Checks whether a password entered by the user is correct."""
        if hmac.compare_digest(st.session_state["password"], st.secrets["password"]):
            st.session_state["password_correct"] = True
            del st.session_state["password"]  # Don't store the password.
        else:
            st.session_state["password_correct"] = False

    # Return True if the password is validated.
    if st.session_state.get("password_correct", False):
        return True

    # Show input for password.
    st.text_input(
        "Password", type="password", on_change=password_entered, key="password"
    )
    if "password_correct" in st.session_state:
        st.error("😕 Password incorrect")
    return False

def check_feats_look_circular(seq):
    """Checks to see if features look circular."""
    features = seq.features
    for feature in features:
        if not isinstance(feature.location, CompoundLocation):
            continue
        # Check whether the compound feature is actually a wraparound
        wraparound_part_index = None
        rightmost_part = None
        last_part_start = None
        for i, part in enumerate(feature.location.parts):
            if rightmost_part != None and part.start < last_part_start:
                return True
            if rightmost_part == None:
                rightmost_part = i
                last_part_start = part.start
    else:
        return False

def run_webapp():
    if not check_password():
        st.stop()  # Do not continue if check_password is not True.


    st._config.set_option(f"theme.base", "light")
    st.set_page_config(
        page_title="EFM Calculator",
        page_icon=ASSET_LOCATION + "/favicon.ico",
        layout="wide",
        initial_sidebar_state="auto",
        menu_items=None,
    )

    st.markdown("""
    <style>
        .block-container {
            padding-top: 0rem;
            padding-bottom: 0rem;
            max-width: 1500px;
        }
    </style>
    """, unsafe_allow_html=True)

    st.markdown("""
    <style>
        #MainMenu, header, footer {visibility: hidden;}

        /* This code gets the first element on the sidebar,
        and overrides its default styling */
        section[data-testid="stSidebar"] div:first-child {
            top: 0;
            height: 100vh;
        }
    </style>
    """,unsafe_allow_html=True)

    st.markdown( # https://discuss.streamlit.io/t/image-and-text-next-to-each-other/7627/17
        """
        <style>
        .container {
            display: flex;
        }
        .logo-text {
            font-weight:700 !important;
            font-size:50px !important;
            color: #000000 !important;
            padding-top: 75px !important;
        }
        .logo-img {
            float:right;
        }
        </style>
        """,
        unsafe_allow_html=True
    )

    st.markdown("""
    <style>
        * {
           overflow-anchor: none !important;
           }
    </style>""", unsafe_allow_html=True)

    collogo,_,colbadge = st.columns([2,1,2], vertical_alignment="bottom")
    with collogo:
        st.markdown(
            f"""
            <div class="container">
                <img class="logo-img" src="data:image/svg+xml;base64,{base64.b64encode(open(ASSET_LOCATION + "/tombstone.svg", "rb").read()).decode()}">
                <p class="logo-text">EFM Calculator</p>
            </div>
            """,
            unsafe_allow_html=True
        )
        st.text(f"Version {version('efmcalculator')}")

    with colbadge:
        st.html(r'<a href="https://github.com/barricklab/efm-calculator2"><img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/barricklab/efm-calculator2?style=social&label=barricklab%2Fefm-calculator2"></a>')


    col1,col2,col3 = st.columns([2,1,2])

    with col1:
        upload_option = "Upload files (FASTA, GenBank, or CSV)"
        enter_option = "Copy/Paste Plain Text"
        example_option = "Example"

        option = st.radio(
            "Choose method of submitting sequence:",
            [example_option, upload_option, enter_option],
        )

        inSeq = None

    with col3:
        st.write("The EFM Calculator predicts mutational hotspots as a result of DNA polymerase slippage. It classifies these hotspots into three categories, Short Sequence Repeats, Short Repeated Sequences, and Repeat Mediated Deletions. For more information, please see the paper. If you have found this tool helpful, please remember to cite it as well.")
        st.write("Jack, B. R., Leonard, S. P., Mishler, D. M., Renda, B. A., Leon, D., Suárez, G. A., & Barrick, J. E. (2015). Predicting the Genetic Stability of Engineered DNA Sequences with the EFM Calculator. ACS Synthetic Biology, 4(8), 939–943. https://doi.org/10.1021/acssynbio.5b00068")


    # Initialize session state
    if not st.session_state.get("statemachine", False):
        st.session_state["statemachine"] = StateMachine()
    statemachine = st.session_state["statemachine"]

    with TemporaryDirectory() as tempdir:
        is_circular = False
        if option == upload_option:
            with col1:
                is_circular = st.checkbox(label="Circular Prediction", value=True)
            upload_disclaimer = f"Total sequence length must be less than {MAX_SIZE+1}. CSV files must have a 'seq' column and may have a 'name' column."
            uploaded_files = st.file_uploader("Choose a file:", type=VALID_EXTS, accept_multiple_files = True)
            st.write(upload_disclaimer)
            if uploaded_files:
                inSeq = []
                for uploaded_file in uploaded_files:
                    original_filename = os.path.splitext(uploaded_file.name)[0]
                    uploaded_filetype = Path(uploaded_file.name).suffix
                    # Hash the file to create a safe name
                    filename = Path(
                        tempdir + str(hash(uploaded_file.name)) + uploaded_filetype
                    )
                    with open(filename, "wb") as f:
                        f.write(uploaded_file.getbuffer())
                    sequences = parse_file(filename, use_filename=False, iscircular = is_circular)
                    file_sequences = []

                    for sequence in sequences:
                        filename = original_filename
                        if not sequence.name:
                            sequence.name = f"{filename}"
                        if not sequence.description or sequence.description == '':
                            sequence.description = f"{filename}"
                        file_sequences.append(sequence)

                    inSeq.extend(file_sequences)

                st.success("Files uploaded.")

        elif option == enter_option:
            with col1:
                is_circular = st.checkbox(label="Circular Prediction", value=True)
            upload_disclaimer = f"""<div>
            <p>Total sequence length must be less than {MAX_SIZE+1}.</p>
            </div>"""
            field = st.text_area("Input sequence here:", max_chars=MAX_SIZE)
            st.markdown(upload_disclaimer, unsafe_allow_html=True)
            field = field.replace("\n", "")
            field = field.replace(" ", "")
            field = "".join([i for i in field if not i.isdigit()])
            if field:
                record = SeqRecord(Seq(field), id="sequence")
                record = EFMSequence(record, is_circular, originhash = hashlib.md5(("string" + field).encode()))
                inSeq = [SeqRecord(Seq(field), id="sequence")]

        elif option == example_option:
            with col1:
                gbs = []
                examples_path = "examples/"
                for infile_loc in glob.glob(os.path.join(examples_path, "*.gb")) + glob.glob(os.path.join(examples_path, "*.fasta")):
                    gbs.append(infile_loc.split("/")[-1])
                exampleFile = st.radio("Choose example file:", gbs)
                filepath = Path(examples_path + f"{exampleFile}")
                if filepath:
                    inSeq = parse_file(filepath)

                is_circular = True

        if not inSeq:
            st.stop()
        else:
            pass

        statemachine.import_sequences(inSeq)

        if len(inSeq) == 1:
            disable_dropdown = True
        else:
            disable_dropdown = False

        col4,col5,col6 = st.columns([2,1,2])

        with col4:
            selected_sequence = st.selectbox(
                "Sample:", statemachine.named_sequences.keys(),
                disabled = disable_dropdown
            )
            selectedhash = statemachine.named_sequences[selected_sequence]
            seq_record = statemachine.user_sequences[selectedhash]

        with col6:
            st.write("\n")
            with TemporaryDirectory() as tempdir:
                #bulk_output(results, inSeq, tempdir, skip_vis = True)

                filestream=io.BytesIO() # https://stackoverflow.com/questions/75304410/streamlit-download-button-not-working-when-trying-to-download-files-as-zip
                with zipfile.ZipFile(filestream, mode='w', compression=zipfile.ZIP_DEFLATED) as zipf:
                    for root, dirs, files in os.walk(tempdir):
                            for file in files:
                                zipf.write(os.path.join(root, file),
                                            os.path.relpath(os.path.join(root, file),
                                                            os.path.join(tempdir, '..')))
                st.download_button(
                    label="Download results (zip)",
                    data=filestream,
                    file_name="results.zip",
                    mime="application/zip", type="primary"
                )


        unique_features = seq_record.unique_annotations

        if not seq_record.predicted:
            with st.spinner("Calculating..."):
                seq_record.call_predictions(strategy="pairwise")

        figcontainer = st.container(height=640)


        if unique_features:
            feature_filter = st.multiselect('Filter by feature annotation',
                                            sorted(unique_features),
                                            default=seq_record._last_filters)
        else:
            feature_filter = []
        seq_record.set_filters(feature_filter)

        results = [seq_record._filtered_ssrs, seq_record._filtered_srss, seq_record._filtered_rmds]

        ssr_columns = results[0].columns
        srs_columns = results[1].columns
        rmd_columns = results[2].columns
        top_columns = seq_record._filtered_top.columns
        del(top_columns[top_columns.index("show")])
        del(ssr_columns[ssr_columns.index("show")])
        del(srs_columns[srs_columns.index("show")])
        del(rmd_columns[rmd_columns.index("show")])

        summary = rip_score(results[0], results[1], results[2], sequence_length = len(seq_record.seq))
        looks_circular = check_feats_look_circular(seq_record)
        if looks_circular:
            st.warning("You deselected the circular option, but your file looks circular.", icon="⚠️")

        tab1, tab2, tab3, tab4 = st.tabs(["Top", "SSR", "SRS", "RMD"])
        with tab1:
            top_table = st.data_editor(seq_record._filtered_top.to_pandas().style.format({"mutation_rate": "{:,.2e}"}),
                                       disabled=top_columns,
                                       hide_index=True,
                                       on_change = seq_record.upate_top_session,
                                       key="topchanges")
        with tab2:
            ssrtable = results[0].to_pandas().style.format({"mutation_rate": "{:,.2e}"})
            st.data_editor(ssrtable,
                          hide_index=True,
                          disabled=ssr_columns,
                          use_container_width=True,
                          key="ssrchanges",
                          on_change = seq_record.update_ssr_session)
        with tab3:
            srstable = results[1].to_pandas().style.format({"mutation_rate": "{:,.2e}"})
            st.data_editor(srstable, hide_index=True, disabled=srs_columns,
            use_container_width=True,
            key="srschanges",
            on_change = seq_record.update_srs_session)
        with tab4:
            rmdtable = results[2].to_pandas().style.format({"mutation_rate": "{:,.2e}"})
            st.data_editor(rmdtable, hide_index=True, disabled=rmd_columns,
            use_container_width=True,
            key="rmdchanges",
            on_change = seq_record.update_rmd_session)

        fig = bokeh_plot(seq_record)
        with figcontainer:
            st.bokeh_chart(fig, use_container_width=True)

    add_vertical_space(4)
