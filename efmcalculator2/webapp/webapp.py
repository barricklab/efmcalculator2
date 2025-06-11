from os.path import normpath
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
from ..constants import VALID_EXTS, MAX_SIZE
from .bokeh_plot import bokeh_plot
from ..pipeline.primary_pipeline import predict
from ..ingest import parse_file, validate_sequences, BadSequenceError

from bokeh.embed import file_html
import streamlit.components.v1 as components
from tempfile import TemporaryDirectory
from pathlib import Path
from streamlit_extras.stylable_container import stylable_container
from streamlit_extras.add_vertical_space import add_vertical_space
from ..pipeline.mutation_rates import rip_score
from .vis_utils import eval_top
from ..StateMachine import StateMachine
from ..ingest import EFMSequence
from time import sleep

from streamlit_javascript import st_javascript

import hashlib

ASSET_LOCATION = os.path.join(os.path.dirname(__file__), "assets")

import pandas as pd
import base64
import json

@st.fragment
def downloadfragment():
    def download_data(tempdir):
        outputdir = tempdir + "/results"
        os.mkdir(outputdir)
        statemachine = st.session_state["statemachine"]
        filetype = st.session_state["dlft"]
        statemachine.save_results(outputdir, filetype=filetype, prediction_style="pairwise")
        with zipfile.ZipFile(f"{outputdir}/results.zip", mode='w', compression=zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(outputdir):
                    for file in files:
                        zipf.write(os.path.join(root, file),
                                    os.path.relpath(os.path.join(root, file),
                                                    os.path.join(outputdir, '..')))
        return f"{outputdir}/results.zip"


    with TemporaryDirectory() as tempdir:
        if not st.session_state.get("dlft"):
            st.session_state["dlft"] = "csv"
        with open(download_data(tempdir), "rb") as fp:
            st.download_button(label = "Download",
                data = fp,
                mime="application/zip",
                use_container_width=True,
                file_name="results.zip")
        st.selectbox("Download File Format",["csv", "parquet"], key="dlft",
            help="CSV files are easily usable in most spreadsheat programs but lack annotation information. Parquet files require specialized tooling but include annotation metadata")


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

    # Initialize session state
    if not st.session_state.get("statemachine", False):
        st.session_state["statemachine"] = StateMachine()
    statemachine = st.session_state["statemachine"]

    col1,col2,col3 = st.columns([2,1,2])
    with col1:
        st.html(r'<a href="https://github.com/barricklab/efmcalculator2"><img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/barricklab/efmcalculator2?style=social&label=barricklab%2Fefmcalculator2"></a>')
        upload_option = "Upload files (FASTA, GenBank, or CSV)"
        enter_option = "Copy/Paste Plain Text"
        example_option = "Example"

        option = st.radio(
            "Choose method of submitting sequence:",
            [example_option, upload_option, enter_option],
        )

        inSeq = None

    with col3:
        try:
            from .._version import version_tuple
            st.markdown(f"Version {version_tuple[0]}.{version_tuple[1]}.{version_tuple[2]} ([{str(version_tuple[4])[1:8]}](https://www.github.com/barricklab/efmcalculator2/commit/{str(version_tuple[4]).split('.')[0][1:8]}))")
        except:
            pass

    with TemporaryDirectory() as tempdir:
        is_circular = True
        field = None
        if option == upload_option:
            with col3:
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
            else:
                st.session_state["statemachine"] = StateMachine()

        elif option == enter_option:
            with col1:
                is_circular = st.checkbox(label="Circular Prediction", value=True)
            with col3:
                upload_disclaimer = f"""<div>
                <p>Total sequence length must be less than {MAX_SIZE+1}.</p>
                </div>"""
                field = st.text_area("Input sequence here:", max_chars=MAX_SIZE)
                st.markdown(upload_disclaimer, unsafe_allow_html=True)
                field = field.replace("\n", "")
                field = field.replace(" ", "")
                field = "".join([i for i in field if not i.isdigit()])
                last_text_input = st.session_state.get("last_text_input", "")
                if field:
                    record = SeqRecord(Seq(field), id="sequence")
                    originhash = hashlib.md5(("string" + field).encode())
                    record = EFMSequence(record, is_circular, originhash)
                    inSeq = [record]
                else:
                    st.session_state["statemachine"] = StateMachine()

        elif option == example_option:
            with col3:
                gbs = []
                examples_path = "examples/"
                for infile_loc in glob.glob(os.path.join(examples_path, "*.gb")) + glob.glob(os.path.join(examples_path, "*.fasta")):
                    gbs.append(infile_loc.split("/")[-1])
                gbs = sorted(gbs)
                exampleFile = st.radio("Choose example file:", gbs)
                filepath = Path(examples_path + f"{exampleFile}")
                if filepath:
                    inSeq = parse_file(filepath, iscircular = True)
                st.write("Examples run in linear mode")


        if not inSeq:
            st.stop()

        for seq in inSeq:
            seq.oneindex = True


        if option == enter_option and field == st.session_state.get("last_text_input", ""):
            pass
        else:
            statemachine.import_sequences(inSeq, max_size=50000, webapp = True)
        st.session_state["last_text_input"] = field

        for seq in statemachine.sequencestates.values():
            if not seq.predicted:
                seq.efmsequence.call_predictions(strategy="pairwise")
                seq.post_predict_processing()
                seq.reset_selected_predictions()

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
            seq_record = statemachine.sequencestates[selectedhash]

        unique_features = seq_record.unique_annotations

        with col6:
            downloadfragment()

        figcontainer = st.container(height=340)

        if unique_features:
            feature_filter = st.multiselect('Filter by feature annotation',
                                            sorted(unique_features),
                                            default=seq_record._last_filters)
        else:
            feature_filter = []
        seq_record.set_filters(feature_filter)

        if not st.session_state.get("last_filter", []) == feature_filter:
            seq_record.reset_selected_predictions()
        st.session_state["last_filter"] = feature_filter

        results = [seq_record._filtered_ssrs, seq_record._filtered_srss, seq_record._filtered_rmds]

        if feature_filter:
            sequence_of_interest = seq_record.annotation_coverage(feature_filter)
        else:
            sequence_of_interest  = len(seq_record.seq)

        summary = rip_score(results[0], results[1], results[2], sequence_length = sequence_of_interest)
        looks_circular = check_feats_look_circular(seq_record)
        if looks_circular:
            st.warning("You deselected the circular option, but your file looks circular.", icon="⚠️")
        col7, col8 = st.columns(2)
        with col7:
            st.markdown(f"<div style='text-align: center;'>RIP score: {summary['rip']:.2f}</div>", unsafe_allow_html=True)
        with col8:
            if summary['ssr_sum'] > 0:
                st.markdown(f"<div style='text-align: center;'>SSRs: {summary['ssr_sum']:.2e}</div>", unsafe_allow_html=True)
            else:
                st.markdown(f"<div style='text-align: center;'>SSRs: 0</div>", unsafe_allow_html=True)
            if summary['srs_sum'] > 0:
                st.markdown(f"<div style='text-align: center;'>SRSs: {summary['srs_sum']:.2e}</div>", unsafe_allow_html=True)
            else:
                st.markdown(f"<div style='text-align: center;'>SRSs: 0</div>", unsafe_allow_html=True)
            if summary['rmd_sum'] > 0:
                st.markdown(f"<div style='text-align: center;'>RMDs: {summary['rmd_sum']:.2e}</div>", unsafe_allow_html=True)
            else:
                st.markdown(f"<div style='text-align: center;'>RMDs: 0</div>", unsafe_allow_html=True)
            if summary['bps_sum'] > 0:
                st.markdown(f"<div style='text-align: center;'>Base Pair Substitution Rate: {summary['bps_sum']:.2e}</div>", unsafe_allow_html=True)
            else:
                st.markdown(f"<div style='text-align: center;'>Base Pair Substitution Rate: 0</div>", unsafe_allow_html=True)

        tab1, tab2, tab3, tab4 = st.tabs(["Top", "SSR", "SRS", "RMD"])
        seq_record.refresh_last_shown()
        seq_record.rebuild_top_table()
        seq_record.rebuild_ssr_table()
        seq_record.rebuild_srs_table()
        seq_record.rebuild_rmd_table()

        with tab1:
            seq_record.top_webapp_table

        with tab2:
            seq_record.ssr_webapp_table

        with tab3:
            seq_record.srs_webapp_table

        with tab4:
            seq_record.rmd_webapp_table

        with figcontainer:
            fig = bokeh_plot(seq_record)
            sleep(0.1) # Helps to avoid two plots displayed
            st.bokeh_chart(fig, use_container_width=True)

    add_vertical_space(4)
    seq_record.refreshed = False

    st.write("The EFM Calculator predicts mutational hotspots as a result of DNA polymerase slippage. It classifies these hotspots into three categories, Short Sequence Repeats, Short Repeated Sequences, and Repeat Mediated Deletions. For more information, please see the paper. If you have found this tool helpful, please remember to cite it as well.")
    st.write("Jack, B. R., Leonard, S. P., Mishler, D. M., Renda, B. A., Leon, D., Suárez, G. A., & Barrick, J. E. (2015). Predicting the Genetic Stability of Engineered DNA Sequences with the EFM Calculator. ACS Synthetic Biology, 4(8), 939–943. https://doi.org/10.1021/acssynbio.5b00068")
