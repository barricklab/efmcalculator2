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
from bokeh.embed import components
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.core.validation import silence
from bokeh.core.validation.warnings import MISSING_RENDERERS
from ..constants import COLORS
from .features import plot_features
from Bio import SeqFeature
from typing import List
import polars as pl
from rich import print
from .table import generate_bokeh_table

import logging


def draw_rmd(fig, rmd_df):
    rmd_y_pos = 1200
    rmd_shape = ((0, 0), (1, 0), (1, 350), (0, 350), (0, 0))
    return _draw_rmd_logic(fig, rmd_df, rmd_shape, "RMD")


def _draw_rmd_logic(fig, df, shape, type):
    rmd_y_pos = 1500
    columns = df.columns
    ssr_source = {
        "x": [],
        "y": [],
        "color": [],
        "name": [],
        "position": [],
        "mutation_rate": [],
    }

    def map_function(row_results):
        ssr_size = row_results[columns.index("repeat_len")]
        start_position = row_results[columns.index("position_left")]
        drawn_ssr = [x[0] * ssr_size + start_position for x in shape]
        ssr_ys = [x[1] + rmd_y_pos for x in shape]
        ssr_source["x"].append(drawn_ssr)
        ssr_source["y"].append(ssr_ys)
        ssr_source["color"].append("black")
        ssr_source["name"].append(type)
        ssr_source["position"].append(start_position)
        ssr_source["mutation_rate"].append(row_results[columns.index("mutation_rate")])

        return 1

    df.map_rows(
        map_function,
    )

    ssr_glyphs = fig.patches(
        "x",
        "y",
        color="color",
        source=ssr_source,
        alpha=0.5,
        line_color="black",
        line_width=1,
    )

    ssr_glyphs_hover = HoverTool(renderers=[ssr_glyphs], tooltips=[("Name", "@name")])

    fig.add_tools(ssr_glyphs_hover)
    table = generate_bokeh_table(df, type)

    return fig, table
