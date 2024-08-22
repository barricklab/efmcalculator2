import streamlit as st
import subprocess
import io
import glob
import Bio.SeqIO as SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import base64
import zipfile

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
from ..efmcalculator import bulk_output, predict_many

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


    st.markdown("""
    <style>
        .block-container {
            padding-top: 0rem;
            padding-bottom: 0rem;
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

    with colbadge:
        st.html(r'<a href="https://github.com/barricklab/efm-calculator2"><img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/barricklab/ostir?style=social&label=barricklab%2FOSTIR"></a>')


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

    with TemporaryDirectory() as tempdir:

        if option == upload_option or option == enter_option:
            with col1:
                license_accepted = st.checkbox("Accept License", value=False, label_visibility="visible")

        else:
            license_accepted = True

        if not license_accepted:
            with col1:
                st.write("""The EFM Calculator is availablue under the [UT Austin Research License](https://github.com/barricklab/efm-calculator2?tab=License-1-ov-file). This tool is free for use in
                            academic, research, experimental and personal use (but specifically excluding Commercial Use). By continuing, you agree to
                            abide by the terms of this license.""")
            return

        if option == upload_option:
            upload_disclaimer = f"Total sequence length must be less than {MAX_SIZE+1}. CSV files must have a 'seq' column and may have a 'name' column."
            uploaded_files = st.file_uploader("Choose a file:", type=VALID_EXTS, accept_multiple_files = True)
            st.write(upload_disclaimer)
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
            field = st.text_area("Input sequence here:", max_chars=MAX_SIZE)
            st.markdown(upload_disclaimer, unsafe_allow_html=True)
            field = field.replace("\n", "")
            field = field.replace(" ", "")
            field = "".join([i for i in field if not i.isdigit()])
            if field:
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

        if inSeq:
            validate_sequences(inSeq, MAX_SIZE)
            with st.spinner("Calculating..."):
                results = predict_many(
                    sequences = inSeq,
                    strategy = "pairwise",
                    isCircular = True,
                )


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

            col4,col5,col6 = st.columns([2,1,2])

            with col4:
                selected_sequence = st.selectbox(
                    "Sample:", sequence_names,
                    disabled = disable_dropdown
                )

            with col6:
                st.write("\n")
                with TemporaryDirectory() as tempdir:
                    bulk_output(results, inSeq, tempdir, skip_vis = True)

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

            seq_record = sequence_dict[selected_sequence]
            sequence = seq_record.seq.strip("\n\n").upper().replace("U", "T")

            with st.spinner("Calculating..."):
                ssr, srs, rmd = predict(sequence, strategy="pairwise", isCircular=False)
                ssr, srs, rmd = post_process(ssr, srs, rmd)
                result = [ssr, srs, rmd]

                fig, tables = make_plot(seq_record, ssr=result[0], srs=result[1], rmd=result[2])
                # layout = make_webpage(fig, tables)
                layout = add_tables(fig, tables)
                shown_result = components.html(file_html(layout, "cdn"), height=650)
