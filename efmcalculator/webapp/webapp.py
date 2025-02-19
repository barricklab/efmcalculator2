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

from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode

import hashlib

ASSET_LOCATION = os.path.join(os.path.dirname(__file__), "../visualization/assets")

import pandas as pd
import base64
import json



def download_data(): # https://gist.github.com/snehankekre/2dcce9fb42b2f7e1742de7431326b263
    with TemporaryDirectory() as tempdir:
        outputdir = tempdir + "/results"
        os.mkdir(outputdir)
        statemachine = st.session_state["statemachine"]
        filetype = st.session_state["dlft"]
        statemachine.save_results(outputdir, filetype=filetype)
        filestream=io.BytesIO() # https://stackoverflow.com/questions/75304410/streamlit-download-button-not-working-when-trying-to-download-files-as-zip
        with zipfile.ZipFile(filestream, mode='w', compression=zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(outputdir):
                    for file in files:
                        zipf.write(os.path.join(root, file),
                                    os.path.relpath(os.path.join(root, file),
                                                    os.path.join(outputdir, '..')))
        b64 = base64.b64encode(filestream.getvalue()).decode()

    dl_link = f"""
    <html>
    <head>
    <title>Start Auto Download file</title>
    <script src="http://code.jquery.com/jquery-3.2.1.min.js"></script>
    <script>
    $('<a href="data:application/zip;base64,{b64}" download="results.zip">')[0].click()
    </script>
    </head>
    </html>
    """
    components.html(
        dl_link,
        height=0,
    )

def check_password():
    """Returns `True` if the user had the correct password."""
    def password_entered():
        """Checks whether a password entered by the user is correct."""
        if hmac.compare_digest(st.session_state.get("password", ""), st.secrets["password"]):
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
        pass

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
        try:
            from .._version import version_tuple
            st.markdown(f"Version {version_tuple[0]}.{version_tuple[1]}.{version_tuple[2]} ([{str(version_tuple[4])[1:8]}](https://www.github.com/barricklab/efm-calculator2/commit/{str(version_tuple[4]).split('.')[0][1:8]}))")
        except:
            pass
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
        is_circular = True
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
                originhash = hashlib.md5(("string" + field).encode())
                record = EFMSequence(record, is_circular, originhash)
                inSeq = [record]

        elif option == example_option:
            with col1:
                gbs = []
                examples_path = "examples/"
                for infile_loc in glob.glob(os.path.join(examples_path, "*.gb")) + glob.glob(os.path.join(examples_path, "*.fasta")):
                    gbs.append(infile_loc.split("/")[-1])
                exampleFile = st.radio("Choose example file:", gbs)
                filepath = Path(examples_path + f"{exampleFile}")
                if filepath:
                    inSeq = parse_file(filepath, iscircular = True)


        if not inSeq:
            st.stop()

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
                with TemporaryDirectory() as tempdir:
                    submit = st.button("Download results",
                                       on_click=download_data,
                                       use_container_width=True,
                                       type="primary")
                st.selectbox("Download File Format",["csv", "parquet"], key="dlft",
                    help="CSV files are easily usable in most spreadsheat programs but lack annotation information. Parquet files require specialized tooling but include annotation metadata")

        unique_features = seq_record.unique_annotations

        if not seq_record.predicted:
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
        ssr_order = ["show"] + ssr_columns
        srs_order = ["show"] + srs_columns
        rmd_order = ["show"] + rmd_columns
        top_order = ["show"] + top_columns

        column_config = {"predid": None,
            "annotationobjects": None,
            "show": st.column_config.CheckboxColumn("Show", help="Show/hide this prediction"),
            "repeat": st.column_config.TextColumn("Sequence", help="One monomer of the repeated sequence"),
            "source": st.column_config.TextColumn("Classification", help="Is the repeat an SSR, SRS, or RMD?"),
            "repeat_len": st.column_config.NumberColumn("Repeat Length", help="The length of a monomer of the repeat"),
            "start": st.column_config.NumberColumn("Start", help="The start position of the repeat"),
            "count": st.column_config.NumberColumn("Count", help="The number of monomers in the repeat"),
            "first_repeat": st.column_config.NumberColumn("First Repeat", help="The first repeat in the sequence"),
            "second_repeat": st.column_config.NumberColumn("Second Repeat", help="The second repeat in the sequence"),
            "distance": st.column_config.NumberColumn("Distance", help="The distance between the repeats"),
            "mutation_rate": st.column_config.NumberColumn(
                "Mutation Rate", format="%.2e"
            ),
            "annotations": st.column_config.ListColumn("Annotations", help="The anotations a deletion might truncate or remove.")}

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
                st.markdown(f"<div style='text-align: center;'>Basal mutation rate: {summary['bps_sum']:.2e}</div>", unsafe_allow_html=True)
            else:
                st.markdown(f"<div style='text-align: center;'>Basal mutation rate: 0</div>", unsafe_allow_html=True)

        tab1, tab2, tab3, tab4 = st.tabs(["Top", "SSR", "SRS", "RMD"])
        with tab1:
            top_table = st.data_editor(seq_record._filtered_top,
                                       disabled=top_columns,
                                       hide_index=True,
                                       on_change = seq_record.upate_top_session,
                                       key="topchanges",
                                       column_config=column_config,
                                       column_order=top_order,
                                       use_container_width=True)
        with tab2:
            ssrtable = results[0]
            #st.data_editor(ssrtable,
            #              hide_index=True,
            #              disabled=ssr_columns,
            #              use_container_width=True,
            #              key="ssrchanges",
            #              on_change = seq_record.update_ssr_session,
            #              column_config=column_config,
            #             column_order=ssr_order)
            pdssrtable = ssrtable.to_pandas()
            #print(ssrtable)
            #print(pdssrtable)
            #pdssrtable["repeat"] = pdssrtable["repeat"].astype(str)
            #pdssrtable["annotations"] = pdssrtable["annotations"].astype(str)
            
            builder = GridOptionsBuilder.from_dataframe(pdssrtable)
            
            builder.configure_column("show", header_name="Show", cellEditor="agCheckboxCellEditor", editable=True)
            builder.configure_column("repeat", header_name="Sequence", tooltipField="repeat")
            builder.configure_column("repeat_len", header_name="Repeat Length", type=["numericColumn"])
            builder.configure_column("start", header_name="Start", type=["numericColumn"])
            builder.configure_column("count", header_name="Count", type=["numericColumn"])
            builder.configure_column("mutation_rate", header_name="Mutation Rate",
                        type=["numericColumn"], valueFormatter="x.toExponential(2)")
            builder.configure_column("annotations", header_name="Annotations", tooltipField="annotations")

            grid_options = builder.build()
            response = AgGrid(pdssrtable, gridOptions=grid_options, height=500, fit_columns_on_grid_load=True, return_edited_values=True,
            data_return_mode='AS_INPUT', update_mode=GridUpdateMode.VALUE_CHANGED, key = "ssrchanges")

            if response:
                updated_df = response["data"]
                if not pdssrtable.equals(updated_df):
                    #converting pandas back to polars to see if it works, a more elegant solution probably exists
                    ssrtable = pl.from_pandas(updated_df)
                    seq_record.update_ssr_session()
                    print("updated ssr table")



        with tab3:
            srstable = results[1]
            #st.data_editor(srstable, hide_index=True, disabled=srs_columns,
            #use_container_width=True,
            #key="srschanges",
            #on_change = seq_record.update_srs_session,
            #column_config=column_config,
            #column_order=srs_order)
            pdsrstable = srstable.to_pandas()
            print(pdsrstable)
            print(pdsrstable)
            pdsrstable["repeat"] = pdsrstable["repeat"].astype(str)
            pdsrstable["annotations"] = pdsrstable["annotations"].astype(str)
            
            builder = GridOptionsBuilder.from_dataframe(pdsrstable)
            
            builder.configure_column("show", header_name="Show", cellEditor="agCheckboxCellEditor", editable=True)
            builder.configure_column("repeat", header_name="Sequence", tooltipField="repeat")
            builder.configure_column("repeat_len", header_name="Repeat Length", type=["numericColumn"])
            builder.configure_column("start", header_name="Start", type=["numericColumn"])
            builder.configure_column("count", header_name="Count", type=["numericColumn"])
            builder.configure_column("mutation_rate", header_name="Mutation Rate",
                        type=["numericColumn"], valueFormatter="x.toExponential(2)")
            builder.configure_column("annotations", header_name="Annotations", tooltipField="annotations")

            grid_options = builder.build()
            AgGrid(pdsrstable, gridOptions=grid_options, height=500, fit_columns_on_grid_load=True)
        
        with tab4:
            rmdtable = results[2]
            #st.data_editor(rmdtable, hide_index=True, disabled=rmd_columns,
            #use_container_width=True,
            #key="rmdchanges",
            #on_change = seq_record.update_rmd_session,
            #column_config=column_config,
            #column_order=rmd_order)
            pdrmdtable = rmdtable.to_pandas()
            print(pdrmdtable)
            pdrmdtable["repeat"] = pdrmdtable["repeat"].astype(str)
            pdrmdtable["annotations"] = pdrmdtable["annotations"].astype(str)
            
            builder = GridOptionsBuilder.from_dataframe(pdrmdtable)
            
            builder.configure_column("show", header_name="Show", cellEditor="agCheckboxCellEditor", editable=True)
            builder.configure_column("repeat", header_name="Sequence", tooltipField="repeat")
            builder.configure_column("repeat_len", header_name="Repeat Length", type=["numericColumn"])
            builder.configure_column("start", header_name="Start", type=["numericColumn"])
            builder.configure_column("count", header_name="Count", type=["numericColumn"])
            builder.configure_column("mutation_rate", header_name="Mutation Rate",
                        type=["numericColumn"], valueFormatter="x.toExponential(2)")
            builder.configure_column("annotations", header_name="Annotations", tooltipField="annotations")

            grid_options = builder.build()
            AgGrid(pdrmdtable, gridOptions=grid_options, height=500, fit_columns_on_grid_load=True)

        with figcontainer:
            fig = bokeh_plot(seq_record)
            sleep(0.1) # Helps to avoid two plots displayed
            st.bokeh_chart(fig, use_container_width=True)

    add_vertical_space(4)
