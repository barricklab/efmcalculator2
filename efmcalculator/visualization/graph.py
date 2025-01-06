import os
import copy
import shutil
import math
import bokeh
import numpy as np
import pandas as pd
from bokeh.transform import linear_cmap
from bokeh.palettes import viridis
from bokeh.layouts import column, row
from bokeh.events import DocumentReady
from bokeh.io import curdoc, export_svg, show
from bokeh.plotting import figure
from bokeh.models import (
    ColumnDataSource,
    Range1d,
    HoverTool,
    LinearAxis,
    TextInput,
    CustomJS,
    Div,
)
from bokeh.model import DataModel
from bokeh.embed import components
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.core.validation import silence
from bokeh.core.validation.warnings import MISSING_RENDERERS
from ..constants import COLORS
from .features import plot_features
from Bio import SeqFeature
from typing import List, Dict
import polars as pl
from .ssr import draw_ssr
from .rmd import draw_rmd
from .srs import draw_srs
from .eval_top import eval_top
from .table import generate_empty_table
from .table import generate_nerfed_bokeh_table
from bokeh.models import Range1d

import logging


class StaggerDatabase(DataModel):
    ssr = []
    srs = []
    rmd = []


def make_plot(seqrecord, **repeat_dataframes):
    ssr_df = repeat_dataframes.get("ssr", None)
    srs_df = repeat_dataframes.get("srs", None)
    rmd_df = repeat_dataframes.get("rmd", None)

    hasAnnotations = False
    if seqrecord.annotations:
        hasAnnotations = True

    curdoc().stagger_database = StaggerDatabase()

    # Set up plot

    fig = figure(width=1250, height=250)

    fig.xaxis.axis_label = "Position"
    fig.yaxis.visible = False
    fig.ygrid.visible = False
    fig.xgrid.visible = False
    fig.toolbar.logo = None

    xmax = len(seqrecord.seq)

    fig.x_range = Range1d(0-len(seqrecord.seq)*.05, len(seqrecord.seq)*1.05)

    ssr_df, srs_df, rmd_df = eval_top(ssr_df, srs_df, rmd_df)

    if hasAnnotations:
        fig = plot_features(seqrecord, fig)

    tables = {}
    if isinstance(ssr_df, pl.DataFrame) and not ssr_df.is_empty():
        fig, table = draw_ssr(fig, ssr_df)
        tables["SSR"] = table
    else:
        tables["SSR"] = generate_empty_table("SSR")

    if isinstance(srs_df, pl.DataFrame) and not srs_df.is_empty():
        fig, table = draw_srs(fig, srs_df)
        tables["SRS"] = table
    else:
        tables["SRS"] = generate_empty_table("SRS")

    if isinstance(rmd_df, pl.DataFrame) and not rmd_df.is_empty():
        fig, table = draw_rmd(fig, rmd_df)
        tables["RMD"] = table
    else:
        tables["RMD"] = generate_empty_table("RMD")

    #recreating with only 2 columns in order to avoid issues while concat, will add position column soon
    if hasAnnotations:
        columns_to_keep = ["repeat", "mutation_rate", "start", "count", "position_left", "position_right", "type"]
    else:
        columns_to_keep = ["repeat", "mutation_rate", "type"]
    ssr_selected = addType(ssr_df, "SSR")
    srs_selected = addType(srs_df, "SRS")
    rmd_selected = addType(rmd_df, "RMD")

    combined_df = pl.concat([ssr_selected, srs_selected, rmd_selected], how="diagonal")
    combined_df = extract_columns(combined_df, columns_to_keep)

    if not combined_df.is_empty():
        top_10_combined = combined_df.sort(by="mutation_rate", descending = True)
        tables["Top"] = generate_nerfed_bokeh_table(top_10_combined)
    else:
        tables["Top"] = generate_empty_table("Top")

    fig.line(
        [0, xmax],
        [1000, 1000],
        line_width=2,
        color="black",
    )


    reordered_tables = {
            "Top": tables["Top"],
            "SSR": tables["SSR"],
            "SRS": tables["SRS"],
            "RMD": tables["RMD"]
        }

    return fig, tables

#cleans up infomation for mutations before presenting specific columns
def extract_columns(df, columns):
    if df is None or df.is_empty():
        return pl.DataFrame({col: [] for col in columns})
    selected_data = {col: df[col] if col in df.columns else pl.lit(None) for col in columns}
    df = df.fill_nan(' ')

    df = pl.DataFrame(selected_data)

    return df

#adds mutation rate type (intended to be used for the "Top" table)
def addType(df, mut_type):
    if mut_type == "SSR":
        df = df.with_columns(pl.lit("SSR").alias("type"))
    elif mut_type == "SRS":
        df = df.with_columns(pl.lit("SRS").alias("type"))
    else:
        df = df.with_columns(pl.lit("RMD").alias("type"))
    return df
