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
    rmd_shape = [[0, 0], [1, 0], [1, 350], [0, 350], [0, 0]]
    return _draw_rmd_logic(fig, rmd_df, rmd_shape, "RMD")


def _draw_rmd_logic(fig, df, shape, type):
    columns = df.columns
    ssr_source = {
        "x": [],
        "y": [],
        "color": [],
        "name": [],
        "position": [],
        "mutation_rate": [],
    }
    empty_rmd_source = copy.deepcopy(ssr_source)

    rmd_glyphs = fig.patches(
        "x",
        "y",
        color="color",
        source=ssr_source,
        alpha=0.5,
        line_color="black",
        line_width=1,
    )

    ssr_glyphs_hover = HoverTool(renderers=[rmd_glyphs], tooltips=[("Name", "@name")])

    fig.add_tools(ssr_glyphs_hover)

    def callback(source_table):
        rmd_y_pos = 1200
        rmd_shape = shape

        javascript = f"""
        var cdata = source_table.data
        var rmd_array = source_table.selected.indices
        var new_glyphs = structuredClone(empty_glyph_source)

        var rmd_shape = {rmd_shape}


        for (let i = 0; i < rmd_array.length; i++) {{

            var len = cdata["repeat_len"][rmd_array[i]]
            var pos_left = cdata["position_left"][rmd_array[i]]
            var mutation_rate = cdata["mutation_rate"][rmd_array[i]]


            new_glyphs['x'].push(rmd_shape.map(x => x[0] * len + pos_left))
            new_glyphs['y'].push(rmd_shape.map(x => x[1] + {rmd_y_pos}))
            new_glyphs['position'].push(pos_left)
            new_glyphs['mutation_rate'].push(mutation_rate)
            new_glyphs['color'].push('black')
            new_glyphs['name'].push('RMD')}}


        glyphs.data_source.data['x'] = new_glyphs['x']
        glyphs.data_source.data['y'] = new_glyphs['y']
        glyphs.data_source.data['position'] = new_glyphs['position']
        glyphs.data_source.data['mutation_rate'] = new_glyphs['mutation_rate']
        glyphs.data_source.data['color'] = new_glyphs['color']
        glyphs.data_source.data['name'] = new_glyphs['name']

        glyphs.data_source.change.emit()
        """

        js_callback = CustomJS(
            args=dict(
                source_table=source_table,
                glyphs=rmd_glyphs,
                empty_glyph_source=empty_rmd_source,
            ),
            code=javascript,
        )

        return js_callback

    table = generate_bokeh_table(df, type, callback=callback)

    return fig, table
