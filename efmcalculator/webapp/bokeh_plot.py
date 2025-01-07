import bokeh as bk
import copy
from Bio import SeqFeature
from bokeh.plotting import figure
from bokeh.models import HoverTool
from typing import List

def bokeh_plot(ssr_df, srs_df, rmd_df, seqobj):
    fig = figure(plot_height=600)
    fig.line(x=[0, len(seqobj.seq)], y=[0, 0], line_color="black", line_width=2)
    fig.xaxis.axis_label = "Position"
    fig.yaxis.visible = False
    fig.ygrid.visible = False
    fig.xgrid.visible = False
    fig.toolbar.logo = None

    fig, depth = plot_annotations(fig, seqobj)


    return fig

def plot_annotations(fig, seqobj):

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

    seqobj.features = assign_feature_levels(seqobj.features)
    xmax = len(seqobj.seq)

    annotation_depth = -750
    lowest_annotation_y = -500
    genbank_dictionary = {
            "x": [],
            "y": [],
            "color": [],
            "name": [],
            "position": [],
            "strand": [],
        }
    for genbank_annotation in seqobj.features:
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
        genbank_dictionary["strand"].append(genbank_annotation.location.strand)


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

    lowest_annotation_y -= 250
    fig.line(
        [0, xmax],
        [lowest_annotation_y, lowest_annotation_y],
        line_width=2,
        color="black",
    )

    return fig, lowest_annotation_y

def plot_ssr(fig, plotted_ssrs):
    columns = plotted_ssrs.columns
    ssr_source = {
        "x": [],
        "y": [],
        "color": [],
        "name": [],
        "position": [],
        "sequence": [],
        "mutation_rate": [],
        "line_width": [],
        "line_color": [],
    }
    ssr_outlines = copy.deepcopy(ssr_source)

    ssr_shape = [[0, 0], [1, 0], [1, 1000], [0, 1000], [0, 0]]
    outline_ssr_shape = [[0, 0], [1, 0], [1, 500], [1, 1000], [0, 1000], [0, 500], [0, 0]]






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

    ssr_outlines_hover = HoverTool(
        renderers=[ssr_outlines],
        tooltips=[
            ("Type", "@name"),
            ("Sequence", "@sequence"),
            ("Position", "@position"),
            ("Mutation Rate", "@mutation_rate"),
        ],
    )
    fig.add_tools(ssr_outlines_hover)

    return fig

def plot_rmd(fig, ssr_df):
    ssr_source = {
        "x": [],
        "y": [],
        "color": [],
        "line_color": [],
        "name": [],
        "position_left": [],
        "position_right": [],
        "mutation_rate": [],
        "line_width": [],
        "sequence": [],
        "distance": [],
    }

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

    rmd_line_source = {
        "x": [],
        "y": [],
        "color": [],
    }

    rmd_lines = fig.multi_line(
        xs="x",
        ys="y",
        source=rmd_line_source,
        alpha=0.5,
        color="color",
        line_width=3,
    )

    rmd_glyphs_hover = HoverTool(
        renderers=[rmd_glyphs_outline],
        tooltips=[
            ("Type", "@name"),
            ("Sequence", "@sequence"),
            ("Position (left)", "@position_left"),
            ("Position (right)", "@position_right"),
            ("Distance", "@distance"),
            ("Mutation Rate", "@mutation_rate"),
        ],
    )

    fig.add_tools(rmd_glyphs_hover)


    return fig
