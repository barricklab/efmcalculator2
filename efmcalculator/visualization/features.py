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
from Bio import SeqFeature
from typing import List

import logging

genbank_dictionary = {}
annotation_positions = []
annotation_names = []

def get_annotation_positions():
    return annotation_positions

def get_annotation_names():
    return annotation_names

def plot_features(seqrecord, fig):
    xmax = len(seqrecord.seq)
    annotation_depth = 0

    # annotation_depth -= 500
    lowest_annotation_y = copy.copy(annotation_depth) - 500
    genbank_dictionary = {
        "x": [],
        "y": [],
        "color": [],
        "name": [],
        "position": [],
        "strand": [],
    }

    seqrecord.features = assign_feature_levels(seqrecord.features)

    for genbank_annotation in seqrecord.features:
        arrow_depth = 100
        # @TODO: SCALE ANNOTATIONS

        scaled_start = genbank_annotation.location.start
        scaled_end = genbank_annotation.location.end
        featurename = genbank_annotation.qualifiers.get("label")

        if (
            genbank_annotation.location.end - genbank_annotation.location.start
            < arrow_depth
        ):
            arrow_depth = (
                genbank_annotation.location.end - genbank_annotation.location.start
            )
        annotation_base_y = (
            -1250.0 * genbank_annotation.qualifiers.get("nest_level", 0)
            + annotation_depth
        )

        if genbank_annotation.location.strand == 1:
            xs = [
                scaled_start,
                scaled_start,
                scaled_end - arrow_depth,
                scaled_end,
                scaled_end - arrow_depth,
            ]
            ys = [
                annotation_base_y + 500,
                annotation_base_y - 500,
                annotation_base_y - 500,
                annotation_base_y,
                annotation_base_y + 500,
            ]
        elif genbank_annotation.location.strand == -1:
            xs = [
                scaled_end,
                scaled_end,
                scaled_start + arrow_depth,
                scaled_start,
                scaled_start + arrow_depth,
            ]
            ys = [
                annotation_base_y - 500,
                annotation_base_y + 500,
                annotation_base_y + 500,
                annotation_base_y,
                annotation_base_y - 500,
            ]
        else:
            xs = [scaled_start, scaled_start, scaled_end, scaled_end]
            ys = [
                annotation_base_y - 500,
                annotation_base_y + 500,
                annotation_base_y + 500,
                annotation_base_y - 500,
            ]

        # Update lowest annotation y
        if annotation_base_y - 500 < lowest_annotation_y:
            lowest_annotation_y = annotation_base_y - 500

        # Define the color
        color = get_feature_color(genbank_annotation)

        # Add the annotation to the dictionary

        genbank_dictionary["x"].append(xs)
        genbank_dictionary["y"].append(ys)
        genbank_dictionary["color"].append(color)
        genbank_dictionary["name"].append(featurename)
        genbank_dictionary["position"].append(
            f"{genbank_annotation.location.start}-{genbank_annotation.location.end}"
        )
        #annotation_positions.append(genbank_dictionary["position"])
        genbank_dictionary["strand"].append(genbank_annotation.location.strand)

    annotation_positions.append(genbank_dictionary["position"])
    annotation_names.append(genbank_dictionary["name"])
    #annotation_positions = dict(zip(genbank_dictionary["name"], genbank_dictionary["position"]))
    genbank_glyphs = fig.patches(
        "x",
        "y",
        color="color",
        source=genbank_dictionary,
        alpha=0.5,
        line_color="black",
        line_width=2,
    )
    genbank_glyphs_hover = HoverTool(
        renderers=[genbank_glyphs], tooltips=[("Name", "@name")]
    )
    fig.add_tools(genbank_glyphs_hover)

    # Draw a line below the annotations

    lowest_annotation_y -= 500
    fig.line(
        [0, xmax],
        [lowest_annotation_y, lowest_annotation_y],
        line_width=2,
        color="black",
    )

    return fig

def get_feature_color(feature: SeqFeature) -> str:
    feature_type = feature.type
    if feature.qualifiers.get("ApEinfo_fwdcolor"):
        return feature.qualifiers.get("ApEinfo_fwdcolor")[0]
    elif feature.qualifiers.get("ApEinfo_revcolor"):
        return feature.qualifiers.get("ApEinfo_revcolor")[0]
    elif feature_type in ["gene"]:
        color = COLORS["cds"]
    elif feature_type in ["promoter"]:
        color = COLORS["promoter"]
    elif feature_type in ["terminator"]:
        color = COLORS["terminator"]
    elif feature_type in ["ori", "rep_origin"]:
        color = COLORS["ori"]
    elif feature_type in ["ncRNA"]:
        color = COLORS["ncrna"]
    elif feature_type in ["primer_bind"]:
        color = COLORS["primer_bind"]
    else:
        color = COLORS["misc"]
    return color


def assign_feature_levels(features: List[SeqFeature]) -> List[SeqFeature]:
    end_positions = dict()
    for feature in features:
        if not end_positions:
            end_positions[str(0)] = feature.location.end
            feature.qualifiers["nest_level"] = 0
            continue
        for level, position in end_positions.items():
            level = int(level)
            if position == "inf":
                continue
            elif int(position) > int(feature.location.start):
                continue
            else:
                end_positions[str(level)] = feature.location.end
                feature.qualifiers["nest_level"] = level
                break
        if (
            not feature.qualifiers.get("nest_level")
            and feature.qualifiers.get("nest_level") != 0
        ):  # Layer 0 is falcey
            lowest_level = max([int(n) for n in end_positions])
            end_positions[str(lowest_level + 1)] = feature.location.end
            feature.qualifiers["nest_level"] = lowest_level + 1
    return features
