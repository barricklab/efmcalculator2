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
from bokeh.io import curdoc, export_svg, show, output_file, save
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
from bokeh.models import Styles
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
from ..mutation_rates import rip_score

import logging

FONT = "Arial"
FONTSIZE = "12pt"


def add_header(layout):
    # Build a DIV above the plot that contains the name of the plot and the total burden
    if False:
        name_div = Div(
            text=f'<h1 style="font-family: {FONT}; font-size: {FONTSIZE}">Sample Name</h1>'
        )
    else:
        name_div = Div()
    burden_div = Div(
        text=f'<h2 style="font-family: {FONT}; font-size: {FONTSIZE}">RIP Score</h2>'
    )

    # layout = column(layout, widgets, styles={'margin': '0 auto', 'align-items': 'center'})
    layout = column(
        name_div,
        burden_div,
        layout,
        styles={"margin": "0 auto", "align-items": "center"},
    )

    return layout


def assemble(layout, summary, tables):
    if not layout:
        layout = Div()

    style = Styles(font = "IBM Plex Sans") # <- Doesnt appear to work?
    summary_layout = column(Div(text=f"Rip Score: {summary["rip"]}", styles=style, width=500),
                            Div(text=f"SSR: {summary["ssr_sum"]}", width=500),
                            Div(text=f"SRS: {summary["srs_sum"]}", width=500),
                            Div(text=f"RMD: {summary["rmd_sum"]}", width=500),
                            Div(text=f"Basal: {summary["bps_sum"]}", width=500),)

    if len(tables) > 1:
        # Make buttons that show the appropriat tables
        buttons = []
        button_label = Div(text="Tables:")
        for i, table_data in enumerate(tables.items()):
            table_name = table_data[0]
            table_value = table_data[1]
            button = bokeh.models.widgets.Button(label=table_name)
            button.js_on_click(
                CustomJS(
                    args=dict(tables=tables, label=table_name),
                    code="""
            tables[label].visible = true;
            for (var key in tables) {
                if (key != label){
                    tables[key].visible = false;
                }
            }
            """,
                )
            )
            buttons.append(button)
            if i == 0:
                table_value.visible = True
            else:
                table_value.visible = False

        buttons = row(
            button_label,
            *buttons,
            styles={"margin": "0 auto", "align-items": "center"},
        )

        table_layout = column(
                    buttons, *tables.values(), styles={"margin": "0 auto", }
                )

        data_layout = row(summary_layout, table_layout, styles={"margin": "0 auto", })

        layout = column(
            layout,
            data_layout,
            styles={"margin": "0 auto", "align-items": "center"},
        )

    return layout


def add_footer(layout):
    return layout


def make_standalone_page(fig, tables, summary):
    layout = fig
    layout.sizing_mode = 'stretch_width'

    layout = add_header(layout)
    if tables:
        layout = assemble(layout, summary, tables)
    layout = add_footer(layout)

    return layout


def export_html(layout, filename):
    script, div = components(layout)
    bokeh_version = bokeh.__version__
    template_path = os.path.join(
        os.path.dirname(__file__), "assets", "html_template.html"
    )
    asset_location = os.path.join(os.path.dirname(__file__), "assets")
    with open(template_path, "r", encoding="utf-8") as f:
        template = f.read()

    if os.path.exists(os.path.join(os.path.dirname(filename), "assets")):
        shutil.rmtree(os.path.join(os.path.dirname(filename), "assets"))
    shutil.copytree(
        os.path.join(os.path.dirname(__file__), "assets"),
        os.path.join(os.path.dirname(filename), "assets"),
    )
    template = template.replace("{{bokehscript}}", script)
    template = template.replace("{{bokehdiv}}", div)
    template = template.replace("{{bokehversion}}", bokeh_version)
    template = template.replace("{{assetlocation}}", asset_location)
    with open(filename, "w", encoding="utf-8") as f:
        f.write(template)
    return filename
