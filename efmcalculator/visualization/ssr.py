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


def draw_ssr(fig, ssr_df):
    ssr_y_pos = 500
    ssr_shape = ((0, 0), (1, 0), (1, 350), (0, 350), (0, 0))

    table = generate_bokeh_table(ssr_df, "SSR")

    columns = ssr_df.columns
    ssr_source = {
        "x": [],
        "y": [],
        "color": [],
        "name": [],
        "position": [],
        "mutation_rate": [],
    }

    def map_function(row_results):
        ssr_size = (
            row_results[columns.index("repeat_len")]
            * row_results[columns.index("count")]
        )
        start_position = row_results[columns.index("start")]
        drawn_ssr = [x[0] * ssr_size + start_position for x in ssr_shape]
        ssr_ys = [x[1] + ssr_y_pos for x in ssr_shape]
        ssr_source["x"].append(drawn_ssr)
        ssr_source["y"].append(ssr_ys)
        ssr_source["color"].append("black")
        ssr_source["name"].append("SSR")
        ssr_source["position"].append(start_position)
        ssr_source["mutation_rate"].append(row_results[columns.index("mutation_rate")])

        return 1

    ssr_df.map_rows(
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

    return fig, table
