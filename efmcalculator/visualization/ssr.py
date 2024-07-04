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
    columns = ssr_df.columns
    ssr_source = {
        "x": [],
        "y": [],
        "color": [],
        "name": [],
        "position": [],
        "mutation_rate": [],
    }
    empty_ssr_source = copy.deepcopy(ssr_source)

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

    def callback(source_table):
        ssr_y_pos = 1000
        ssr_shape = [[0, 0], [1, 0], [1, 350], [0, 350], [0, 0]]

        javascript = f"""
        var cdata = source_table.data
        var ssr_array = source_table.selected.indices
        var new_glyphs = structuredClone(empty_glyph_source)

        var ssr_shape = {ssr_shape}


        for (let i = 0; i < ssr_array.length; i++) {{

            var len = cdata["repeat_len"][ssr_array[i]]
            var count = cdata["count"][ssr_array[i]]
            var start = cdata["start"][ssr_array[i]]
            var mutation_rate = cdata["mutation_rate"][ssr_array[i]]


            new_glyphs['x'].push(ssr_shape.map(x => x[0] * len * count + start))
            new_glyphs['y'].push(ssr_shape.map(x => x[1] + {ssr_y_pos}))
            new_glyphs['position'].push(start)
            new_glyphs['mutation_rate'].push(mutation_rate)
            new_glyphs['color'].push('black')
            new_glyphs['name'].push('SSR')}}

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
                glyphs=ssr_glyphs,
                empty_glyph_source=empty_ssr_source,
            ),
            code=javascript,
        )

        return js_callback

    table = generate_bokeh_table(ssr_df, "SSR", callback=callback)

    return fig, table
