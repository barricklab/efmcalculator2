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
        "line_width": [],
        "line_color": [],
    }
    empty_ssr_source = copy.deepcopy(ssr_source)

    ssr_glyphs = fig.patches(
        "x",
        "y",
        color="color",
        source=ssr_source,
        alpha=0.5,
        line_color="black",
        line_width="line_width",
    )
    ssr_outlines = fig.patches(
        "x",
        "y",
        source=ssr_source,
        alpha=0.5,
        color="color",
        line_color="line_color",
        line_width="line_width",
    )

    rmd_glyphs_outline = fig.patches(
        "x",
        "y",
        color="color",
        source=ssr_source,
        alpha=0.5,
        line_color="black",
        line_width=1,
    )

    # ssr_glyphs_hover = HoverTool(renderers=[ssr_glyphs], tooltips=[("Name", "@name")])
    ssr_outlines_hover = HoverTool(
        renderers=[ssr_outlines], tooltips=[("Name", "@name")]
    )
    # fig.add_tools(ssr_glyphs_hover)
    fig.add_tools(ssr_outlines_hover)

    def callback(source_table):
        outline_girth_x = 50
        outline_girth_y = 250
        # [layer_checker, start, (len*count)+start+({outline_girth_x}*2)]

        javascript = f"""
        var cdata = source_table.data
        var ssr_array = source_table.selected.indices
        var new_glyphs = structuredClone(empty_glyph_source)
        var new_outlines = structuredClone(empty_glyph_source)
        var ssr_shape = [[0, 0], [1, 0], [1, 1000], [0, 1000], [0, 0]]
        var outline_ssr_shape = [[0, 0], [1, 0], [1, 500], [1, 1000], [0, 1000], [0, 500], [0, 0]]

        stagger_database.ssr = []
        var mutation_types = ['ssr', 'srs', 'rmd']

        for (let i = 0; i < ssr_array.length; i++) {{
            //Extract information from dataframe
            var len = cdata["repeat_len"][ssr_array[i]]
            var count = cdata["count"][ssr_array[i]]
            var start = cdata["start"][ssr_array[i]]
            var mutation_rate = cdata["mutation_rate"][ssr_array[i]]
            var end = (len*count)+start+{outline_girth_x}

            //Find draw layer
            var layer = -1
            var layer_checker = 1
            var recheck = true
            while (recheck){{
                recheck = false
                for (const property in mutation_types) {{
                    var prop = mutation_types[property]
                    for (var j in stagger_database[prop]) {{
                        if (stagger_database[prop][j][0] != layer_checker){{
                            continue
                        }}

                        if (( stagger_database[prop][j][1]  <=  start-{outline_girth_x} && start-{outline_girth_x} <=  stagger_database[prop][j][2] )
                            || //If start is inside a glyph OR
                            ( stagger_database[prop][j][1]  <=  end && end   <=  stagger_database[prop][j][2] )
                            || // If end is inside a glyph OR
                            (start-{outline_girth_x} <=  stagger_database[prop][j][1] && stagger_database[prop][j][2]  <=  end)){{
                            // If start and end wrap another glyph, then this is an overlapping glyph
                            layer_checker++
                            recheck = true
                            break
                            }}
                    if (recheck){{
                        break
                        }}
                    }}
                }}
            }}
            layer = layer_checker
            stagger_database['ssr'].push([layer, start-{outline_girth_x}, end])
            var ssr_y_pos = 1500*layer

            //Draw glyphs
            var ssr_x = ssr_shape.map(x => x[0] * len * count + start)
            var ssr_y = ssr_shape.map(x => x[1] + ssr_y_pos)
            var ssr_x_outline = outline_ssr_shape.map(x => x[0] * len * count + start)
            var ssr_y_outline = outline_ssr_shape.map(x => x[1] + ssr_y_pos)

            ssr_x_outline[0] = ssr_x_outline[0]-{outline_girth_x}/2
            ssr_x_outline[1] = ssr_x_outline[1]+{outline_girth_x}/2
            ssr_x_outline[2] = ssr_x_outline[2]+{outline_girth_x}
            ssr_x_outline[3] = ssr_x_outline[3]+{outline_girth_x}/2
            ssr_x_outline[4] = ssr_x_outline[4]-{outline_girth_x}/2
            ssr_x_outline[5] = ssr_x_outline[5]-{outline_girth_x}
            ssr_x_outline[6] = ssr_x_outline[6]-{outline_girth_x}/2

            ssr_y_outline[0] = ssr_y_outline[0]-{outline_girth_y}
            ssr_y_outline[1] = ssr_y_outline[1]-{outline_girth_y}
            ssr_x_outline[2] = ssr_x_outline[2]
            ssr_y_outline[3]= ssr_y_outline[3]+{outline_girth_y}
            ssr_y_outline[4] = ssr_y_outline[4]+{outline_girth_y}
            ssr_x_outline[5] = ssr_x_outline[5]
            ssr_y_outline[6] = ssr_y_outline[6]-{outline_girth_y}

            new_glyphs['x'].push(ssr_x)
            new_glyphs['y'].push(ssr_y)
            new_glyphs['position'].push(start)
            new_glyphs['mutation_rate'].push(mutation_rate)
            new_glyphs['color'].push('black')
            new_glyphs['line_color'].push('black')
            new_glyphs['name'].push('SSR')
            new_glyphs['line_width'].push(1)

            new_outlines['x'].push(ssr_x_outline)
            new_outlines['y'].push(ssr_y_outline)
            new_outlines['position'].push(start)
            new_outlines['mutation_rate'].push(mutation_rate)
            new_outlines['color'].push('none')
            new_outlines['line_color'].push('red')
            new_outlines['name'].push('SSR')
            new_outlines['line_width'].push(3)}}

        //Push glyphs to dataframe
        glyphs.data_source.data['x'] = new_glyphs['x']
        glyphs.data_source.data['y'] = new_glyphs['y']
        glyphs.data_source.data['position'] = new_glyphs['position']
        glyphs.data_source.data['mutation_rate'] = new_glyphs['mutation_rate']
        glyphs.data_source.data['color'] = new_glyphs['color']
        glyphs.data_source.data['line_color'] = new_glyphs['line_color']
        glyphs.data_source.data['name'] = new_glyphs['name']
        glyphs.data_source.data['line_width'] = new_glyphs['line_width']


        outlines.data_source.data['x'] = new_outlines['x']
        outlines.data_source.data['y'] = new_outlines['y']
        outlines.data_source.data['position'] = new_outlines['position']
        outlines.data_source.data['mutation_rate'] = new_outlines['mutation_rate']
        outlines.data_source.data['color'] = new_outlines['color']
        outlines.data_source.data['line_color'] = new_outlines['line_color']
        outlines.data_source.data['name'] = new_outlines['name']
        outlines.data_source.data['line_width'] = new_outlines['line_width']

        glyphs.data_source.change.emit()
        outlines.data_source.change.emit()
        """

        js_callback = CustomJS(
            args=dict(
                source_table=source_table,
                glyphs=ssr_glyphs,
                outlines=ssr_outlines,
                empty_glyph_source=empty_ssr_source,
                stagger_database=curdoc().stagger_database,
            ),
            code=javascript,
        )

        curdoc().js_on_event(DocumentReady, js_callback)

        return js_callback

    table = generate_bokeh_table(ssr_df, "SSR", callback=callback)

    return fig, table
