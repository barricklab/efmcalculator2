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
    rmd_color = "yellow"
    return _draw_rmd_logic(fig, rmd_df, rmd_color, "RMD")


def _draw_rmd_logic(fig, df, color, type):
    columns = df.columns
    ssr_source = {
        "x": [],
        "y": [],
        "color": [],
        "line_color": [],
        "name": [],
        "position": [],
        "mutation_rate": [],
        "line_width": [],
    }
    empty_rmd_source = copy.deepcopy(ssr_source)

    rmd_glyphs_left = fig.patches(
        "x",
        "y",
        color="color",
        source=ssr_source,
        alpha=0.5,
        line_color="line_color",
        line_width="line_width",
    )

    rmd_glyphs_right = fig.patches(
        "x",
        "y",
        color="color",
        source=ssr_source,
        alpha=0.5,
        line_color="line_color",
        line_width="line_width",
    )

    rmd_glyphs_outline = fig.patches(
        "x",
        "y",
        color="color",
        source=ssr_source,
        alpha=0.5,
        line_color="line_color",
        line_width="line_width",
    )

    rmd_glyphs_hover = HoverTool(
        renderers=[rmd_glyphs_outline], tooltips=[("Name", "@name")]
    )

    fig.add_tools(rmd_glyphs_hover)

    def callback(source_table):
        rmd_y_pos = 3000

        outline_girth_x = 50
        outline_girth_y = 250

        javascript = f"""
        var cdata = source_table.data
        var rmd_array = source_table.selected.indices
        var new_glyphs_left = structuredClone(empty_glyph_source)
        var new_glyphs_right = structuredClone(empty_glyph_source)
        var new_outlines = structuredClone(empty_glyph_source)

        var rmd_shape = [[0, 0], [1, 0], [1, 1000], [0, 1000], [0, 0]]
        var outline_rmd_shape_squared = structuredClone(rmd_shape)


        for (let i = 0; i < rmd_array.length; i++) {{


            var len = cdata["repeat_len"][rmd_array[i]]
            var pos_left = cdata["position_left"][rmd_array[i]]
            var pos_right = cdata["position_right"][rmd_array[i]]
            var mutation_rate = cdata["mutation_rate"][rmd_array[i]]

            var x_left = rmd_shape.map(x => x[0] * len + pos_left)
            var y = rmd_shape.map(x => x[1] + {rmd_y_pos})
            var x_right = rmd_shape.map(x => x[0] * len + pos_right)

            new_glyphs_left['x'].push(x_left)
            new_glyphs_left['y'].push(y)
            new_glyphs_left['position'].push(pos_left)
            new_glyphs_left['mutation_rate'].push(mutation_rate)
            new_glyphs_left['color'].push('black')
            new_glyphs_left['name'].push('RMD')
            new_glyphs_left['line_color'].push('black')
            new_glyphs_left['line_width'].push(1)

            new_glyphs_right['x'].push(x_right)
            new_glyphs_right['y'].push(y)
            new_glyphs_right['position'].push(pos_right)
            new_glyphs_right['mutation_rate'].push(mutation_rate)
            new_glyphs_right['color'].push('black')
            new_glyphs_right['name'].push('RMD')
            new_glyphs_right['line_color'].push('black')
            new_glyphs_right['line_width'].push(1)

            var x_outline = outline_rmd_shape_squared.map(x => x[0] * len)
            x_outline[0] = x_outline[0]-{outline_girth_x}+pos_left
            x_outline[1] = x_outline[1]+{outline_girth_x}+pos_right
            x_outline[2] = x_outline[2]+{outline_girth_x}+pos_right
            x_outline[3] = x_outline[3]-{outline_girth_x}+pos_left
            x_outline[4] = x_outline[4]-{outline_girth_x}+pos_left

            var y_outline = outline_rmd_shape_squared.map(x => x[1] + {rmd_y_pos})
            y_outline[0] = y_outline[0]-{outline_girth_y}
            y_outline[1] = y_outline[1]-{outline_girth_y}
            y_outline[2] = y_outline[2]+{outline_girth_y}
            y_outline[3] = y_outline[3]+{outline_girth_y}
            y_outline[4] = y_outline[4]-{outline_girth_y}

            new_outlines['x'].push(x_outline)
            new_outlines['y'].push(y_outline)
            new_outlines['position'].push(pos_left)
            new_outlines['mutation_rate'].push(mutation_rate)
            new_outlines['color'].push('none')
            new_outlines['line_color'].push('{color}')
            new_outlines['name'].push('RMD')
            new_outlines['line_width'].push(3)}}



        glyphs_left.data_source.data['x'] = new_glyphs_left['x']
        glyphs_left.data_source.data['y'] = new_glyphs_left['y']
        glyphs_left.data_source.data['position'] = new_glyphs_left['position']
        glyphs_left.data_source.data['mutation_rate'] = new_glyphs_left['mutation_rate']
        glyphs_left.data_source.data['color'] = new_glyphs_left['color']
        glyphs_left.data_source.data['name'] = new_glyphs_left['name']
        glyphs_left.data_source.data['line_color'] = new_glyphs_left['line_color']
        glyphs_left.data_source.data['line_width'] = new_glyphs_left['line_width']

        glyphs_right.data_source.data['x'] = new_glyphs_right['x']
        glyphs_right.data_source.data['y'] = new_glyphs_right['y']
        glyphs_right.data_source.data['position'] = new_glyphs_right['position']
        glyphs_right.data_source.data['mutation_rate'] = new_glyphs_right['mutation_rate']
        glyphs_right.data_source.data['color'] = new_glyphs_right['color']
        glyphs_right.data_source.data['line_color'] = new_glyphs_right['line_color']
        glyphs_right.data_source.data['name'] = new_glyphs_right['name']
        glyphs_right.data_source.data['line_width'] = new_glyphs_right['line_width']

        outlines.data_source.data['x'] = new_outlines['x']
        outlines.data_source.data['y'] = new_outlines['y']
        outlines.data_source.data['position'] = new_outlines['position']
        outlines.data_source.data['mutation_rate'] = new_outlines['mutation_rate']
        outlines.data_source.data['color'] = new_outlines['color']
        outlines.data_source.data['line_color'] = new_outlines['line_color']
        outlines.data_source.data['name'] = new_outlines['name']
        outlines.data_source.data['line_width'] = new_outlines['line_width']

        glyphs_left.data_source.change.emit()
        glyphs_right.data_source.change.emit()
        outlines.data_source.change.emit()
        """

        js_callback = CustomJS(
            args=dict(
                source_table=source_table,
                glyphs_left=rmd_glyphs_left,
                glyphs_right=rmd_glyphs_right,
                outlines=rmd_glyphs_outline,
                empty_glyph_source=empty_rmd_source,
            ),
            code=javascript,
        )

        return js_callback

    table = generate_bokeh_table(df, type, callback=callback)

    return fig, table
