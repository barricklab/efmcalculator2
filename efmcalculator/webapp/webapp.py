import streamlit as st
import subprocess
import io
import glob
import Bio.SeqIO as SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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

    def fix_mobile_columns():  # https://github.com/streamlit/streamlit/issues/6592#issuecomment-1887155722
        st.write(
            """<style>
        [data-testid="column"] {
            width: calc(16.6666% - 1rem) !important;
            flex: 1 1 calc(16.6666% - 1rem) !important;
            min-width: calc(16.6666% - 1rem) !important;
        }
        </style>""",
            unsafe_allow_html=True,
        )

    fix_mobile_columns()

    col1, col2 = st.columns([20, 80])
    with col1:
        st.image(ASSET_LOCATION + "/tombstone.svg", width=150)
    with col2:
        st.title("EFM Calculator")

        upload_option = "Upload a sequence file (FASTA or GenBank)"
        bulk_option = "Upload a CSV"
        enter_option = "Copy/Paste Plain Text"
        example_option = "Example"

        option = st.radio(
            "Choose method of submitting sequence:",
            [upload_option, bulk_option, enter_option, example_option],
        )

    inSeq = None

    with TemporaryDirectory() as tempdir:
        if option == upload_option:
            upload_disclaimer = f"""<div>
            <p>Total sequence length must be less than {MAX_SIZE+1}.</p>
            </div>"""
            st.markdown(upload_disclaimer, unsafe_allow_html=True)
            uploaded_file = st.file_uploader("Choose a file:", type=VALID_EXTS)
            if uploaded_file is not None:
                uploaded_filetype = Path(uploaded_file.name).suffix
                # Hash the file to create a safe name
                filename = Path(
                    tempdir + str(hash(uploaded_file.name)) + uploaded_filetype
                )
                with open(filename, "wb") as f:
                    f.write(uploaded_file.getbuffer())
                sequences = parse_file(filename)
                validate_sequences(sequences, MAX_SIZE)
                st.success("File uploaded.")
                inSeq = sequences[0]

        elif option == bulk_option:
            upload_disclaimer = f"""<div>
            <p>CSV file must have a 'seq' column and may have a 'name' column. Total sequence length must be less than {MAX_SIZE+1}.</p>
            </div>"""
            st.markdown(upload_disclaimer, unsafe_allow_html=True)
            uploaded_file = st.file_uploader("Choose a file:", type="csv")

            if uploaded_file is not None:
                uploaded_filetype = Path(uploaded_file.name).suffix
                # Hash the file to create a safe name
                filename = Path(
                    tempdir + str(hash(uploaded_file.name)) + uploaded_filetype
                )
                with open(filename, "wb") as f:
                    f.write(uploaded_file.getbuffer())
                sequences = parse_file(filename)
                validate_sequences(sequences, MAX_SIZE)
                st.success("File uploaded.")
                inSeq = sequences[0]

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
                inSeq = SeqRecord(Seq(field), id="sequence")
                validate_sequences([inSeq], MAX_SIZE)

        elif option == example_option:
            gbs = []
            examples_path = "examples/"
            for infile_loc in glob.glob(os.path.join(examples_path, "*.gb")):
                gbs.append(infile_loc.split("/")[-1].split(".gb")[0])
            exampleFile = st.radio("Choose example file:", gbs)
            filepath = Path(examples_path + f"{exampleFile}.gb")
            inSeq = parse_file(filepath)
            validate_sequences(inSeq, MAX_SIZE)
            inSeq = inSeq[0]

        if inSeq:
            with st.spinner("Calculating..."):
                sequence = inSeq.seq.strip("\n").upper().replace("U", "T")

                ssr, srs, rmd = predict(sequence, strategy="pairwise", isCircular=False)
                ssr, srs, rmd = post_process(ssr, srs, rmd)
                result = [ssr, srs, rmd]

                fig, tables = make_plot(inSeq, ssr=result[0], srs=result[1], rmd=result[2])
                # layout = make_webpage(fig, tables)
                layout = add_tables(fig, tables)
                components.html(file_html(layout, "cdn"), height=650)

        footer_html = """<div style='text-align: center;'>
        <p>If you've found this tool useful, please cite the following publication: Jack, B. R., Leonard, S. P., Mishler, D. M., Renda, B. A., Leon, D., Suárez, G. A., & Barrick, J. E. (2015). Predicting the Genetic Stability of Engineered DNA Sequences with the EFM Calculator. ACS Synthetic Biology, 4(8), 939–943. https://doi.org/10.1021/acssynbio.5b00068</p>
        </p>
        </div>"""
        st.markdown(footer_html, unsafe_allow_html=True)
