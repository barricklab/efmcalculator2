import streamlit as st
import subprocess
import io
import glob
import Bio.SeqIO as SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import base64

import os
from ..short_seq_finder import predict
from ..visualization.graph import make_plot
from ..visualization.make_webpage import add_tables
from ..constants import VALID_EXTS
from ..parse_inputs import parse_file, validate_sequences, BadSequenceError

from ..efmcalculator import post_process
from bokeh.embed import file_html
import streamlit.components.v1 as components
from tempfile import TemporaryDirectory
from pathlib import Path
from streamlit_extras.stylable_container import stylable_container

ASSET_LOCATION = os.path.join(os.path.dirname(__file__), "../visualization/assets")
MAX_SIZE = 50000


def run_webapp():
    st._config.set_option(f"theme.base", "light")
    st.set_page_config(
        page_title="EFM Calculator",
        page_icon=ASSET_LOCATION + "/favicon.ico",
        layout="wide",
        initial_sidebar_state="auto",
        menu_items=None,
    )

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

    st.markdown(
        f"""
        <div class="container">
            <img class="logo-img" src="data:image/svg+xml;base64,{base64.b64encode(open(ASSET_LOCATION + "/tombstone.svg", "rb").read()).decode()}">
            <p class="logo-text">EFM Calculator</p>
        </div>
        """,
        unsafe_allow_html=True
    )

    upload_option = "Upload files (FASTA, GenBank, or CSV)"
    enter_option = "Copy/Paste Plain Text"
    example_option = "Example"

    option = st.radio(
        "Choose method of submitting sequence:",
        [upload_option, enter_option, example_option],
    )

    inSeq = None

    with TemporaryDirectory() as tempdir:
        if option == upload_option:
            upload_disclaimer = f"Total sequence length must be less than {MAX_SIZE+1}. CSV files must have a 'seq' column and may have a 'name' column."
            st.write(upload_disclaimer)
            uploaded_files = st.file_uploader("Choose a file:", type=VALID_EXTS, accept_multiple_files = True)
            if uploaded_files:
                inSeq = []
                for uploaded_file in uploaded_files:
                    uploaded_filetype = Path(uploaded_files.name).suffix
                    # Hash the file to create a safe name
                    filename = Path(
                        tempdir + str(hash(uploaded_files.name)) + uploaded_filetype
                    )
                    with open(filename, "wb") as f:
                        f.write(uploaded_files.getbuffer())
                    inSeq.extend(parse_file(filename))
                st.success("Files uploaded.")

        elif option == enter_option:
            upload_disclaimer = f"""<div>
            <p>Total sequence length must be less than {MAX_SIZE+1}.</p>
            </div>"""
            st.markdown(upload_disclaimer, unsafe_allow_html=True)
            field = st.text_area("Input sequence here:", max_chars=MAX_SIZE)
            field = field.replace("\n", "")
            field = field.replace(" ", "")
            field = "".join([i for i in field if not i.isdigit()])
            if field:
                inSeq = [SeqRecord(Seq(field), id="sequence")]

        elif option == example_option:
            gbs = []
            examples_path = "examples/"
            for infile_loc in glob.glob(os.path.join(examples_path, "*.gb")) + glob.glob(os.path.join(examples_path, "*.fasta")):
                gbs.append(infile_loc.split("/")[-1])
            exampleFile = st.radio("Choose example file:", gbs)
            filepath = Path(examples_path + f"{exampleFile}")
            if filepath:
                inSeq = parse_file(filepath)


        if inSeq:
            validate_sequences(inSeq, MAX_SIZE)
            with st.spinner("Calculating..."):
                sequence_dict = {}
                sequence_names = []
                for i, seq in enumerate(inSeq):
                    if seq.description:
                        sequence_name = f"{i+1}_{seq.description}"
                    else:
                        sequence_name = f"{i+1}_Sequence"
                    sequence_names.append(sequence_name)
                    sequence_dict[sequence_name] = seq

                if len(inSeq) == 1:
                    disable_dropdown = True
                else:
                    disable_dropdown = False

                selected_sequence = st.selectbox(
                    "Sample:", sequence_names,
                    disabled = disable_dropdown
                )

                seq_record = sequence_dict[selected_sequence]
                sequence = seq_record.seq.strip("\n").upper().replace("U", "T")

                ssr, srs, rmd = predict(sequence, strategy="pairwise", isCircular=False)
                ssr, srs, rmd = post_process(ssr, srs, rmd)
                result = [ssr, srs, rmd]

                fig, tables = make_plot(seq_record, ssr=result[0], srs=result[1], rmd=result[2])
                # layout = make_webpage(fig, tables)
                layout = add_tables(fig, tables)
                shown_result = components.html(file_html(layout, "cdn"), height=650)



        footer_html = """<div style='text-align: center;'>
        <p>If you've found this tool useful, please cite the following publication: Jack, B. R., Leonard, S. P., Mishler, D. M., Renda, B. A., Leon, D., Suárez, G. A., & Barrick, J. E. (2015). Predicting the Genetic Stability of Engineered DNA Sequences with the EFM Calculator. ACS Synthetic Biology, 4(8), 939–943. https://doi.org/10.1021/acssynbio.5b00068</p>
        </p>
        </div>"""
        st.markdown(footer_html, unsafe_allow_html=True)
