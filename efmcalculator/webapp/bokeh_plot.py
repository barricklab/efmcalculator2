from subprocess import run
import bokeh as bk
import copy
from Bio import SeqFeature
from bokeh.plotting import figure
from bokeh.models import HoverTool, CustomJS
from typing import List
import polars as pl
from ..constants import MARKER_HEIGHT, COLORS
from copy import deepcopy

OUTLINE_PADDING_X = 25
OUTLINE_PADDING_Y = 25
INNER_HEIGHT = MARKER_HEIGHT-OUTLINE_PADDING_Y*2
ANNOTATION_HEIGHT = MARKER_HEIGHT/2

def bokeh_plot(seqobj):
    # Figure boilerplate

    fig = figure(plot_height=300, output_backend="svg")
    fig.line(x=[0, len(seqobj.seq)], y=[0, 0], line_color="black", line_width=2)
    fig.xaxis.axis_label = "Position"
    fig.yaxis.visible = False
    fig.ygrid.visible = False
    fig.xgrid.visible = False
    fig.toolbar.logo = None
    fig, depth = plot_annotations(fig, seqobj)

    # Calculate vertical stagger to prevent overlapping annotations
    selected_ssrs = seqobj.get_shown_predictions(seqobj._filtered_ssrs)
    selected_ssrs = seqobj._filtered_ssrs.with_row_index(
                           ).filter(pl.col("index").is_in(selected_ssrs))
    selected_srss = seqobj.get_shown_predictions(seqobj._filtered_srss)
    selected_srss = seqobj._filtered_srss.with_row_index(
                           ).filter(pl.col("index").is_in(selected_srss))
    selected_rmds = seqobj.get_shown_predictions(seqobj._filtered_rmds)
    selected_rmds = seqobj._filtered_rmds.with_row_index(
                           ).filter(pl.col("index").is_in(selected_rmds))

    stagger_database = {
        0: []
    }

    stagger_ssrs = selected_ssrs.with_columns(
        (pl.col("repeat_len")*pl.col("count")+pl.col("start")-1).alias("end").cast(pl.Int32),
        pl.lit("SSR").alias("type"),
        pl.col("start").alias("first_repeat"),
        pl.lit(None).cast(pl.Int32).alias("second_repeat"),
    )
    stagger_srss = selected_srss.with_columns(
        (pl.col("second_repeat") + pl.col("repeat_len")).alias("end").cast(pl.Int32),
        pl.col("first_repeat").alias("start"),
        pl.lit("SRS").alias("type"),
    )
    stagger_rmds = selected_rmds.with_columns(
        (pl.col("second_repeat") + pl.col("repeat_len")).alias("end").cast(pl.Int32),
        pl.col("first_repeat").alias("start"),
        pl.lit("RMD").alias("type"),
    )

    joint_table = pl.concat([stagger_ssrs,
                             stagger_srss,
                             stagger_rmds],
                             how="diagonal").with_columns(
        (pl.col('end') - pl.col('start') + pl.lit(1)).alias("size")
                             ).select(["start", "end", "size", "predid", "type"]).with_columns(
        pl.lit(0).alias("level")
                             )
    joint_table = joint_table.sort(by=["size"], descending=True)

    highest_level = 0
    for i, row in enumerate(joint_table.rows(named=True)):
        assigned = False
        for level in stagger_database:
            for occupied_area in stagger_database[level]:
                if occupied_area[0] <= row['start'] and occupied_area[1] >= row['end']:
                    # Entirely inside
                    break
                elif occupied_area[0] <= row['start'] and row['start'] <= occupied_area[1] <= row['end']:
                    # Left edge
                    break
                elif  row['start'] <= occupied_area[0] <= row['end'] and occupied_area[1] >= row['end']:
                    # Right edge
                    break
            else:
                # Outside all areas
                stagger_database[level].append([row['start']-OUTLINE_PADDING_X*2, row['end']+OUTLINE_PADDING_X*2])
                joint_table[i, 'level'] = level
                assigned = True
                break
        else:
            # All existing levels are full
            highest_level += 1
            stagger_database[highest_level] = [[row['start']-OUTLINE_PADDING_X*2, row['end']+OUTLINE_PADDING_X*2]]
            joint_table[i, 'level'] = highest_level
            assigned = True

    joint_table = joint_table.select(["predid", "level"])
    stagger_ssrs = stagger_ssrs.join(joint_table, on=["predid"], how="left")
    stagger_srss = stagger_srss.join(joint_table, on=["predid"], how="left")
    stagger_rmds = stagger_rmds.join(joint_table, on=["predid"], how="left")

    fig = plot_rmd(fig, stagger_rmds)
    fig = plot_srs(fig, stagger_srss)
    fig = plot_ssr(fig, stagger_ssrs)


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

    staggered_features = assign_feature_levels(seqobj.features)
    xmax = len(seqobj.seq)

    annotation_depth = -500
    lowest_annotation_y = -ANNOTATION_HEIGHT
    genbank_dictionary = {
            "x": [],
            "y": [],
            "color": [],
            "name": [],
            "position": [],
            "strand": [],
        }
    for genbank_annotation in staggered_features:
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
            -650.0 * genbank_annotation.qualifiers.get("nest_level", 0)
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
                annotation_base_y + ANNOTATION_HEIGHT,
                annotation_base_y - ANNOTATION_HEIGHT,
                annotation_base_y - ANNOTATION_HEIGHT,
                annotation_base_y,
                annotation_base_y + ANNOTATION_HEIGHT,
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
                annotation_base_y - ANNOTATION_HEIGHT,
                annotation_base_y + ANNOTATION_HEIGHT,
                annotation_base_y + ANNOTATION_HEIGHT,
                annotation_base_y,
                annotation_base_y - ANNOTATION_HEIGHT,
            ]
        else:
            xs = [scaled_start, scaled_start, scaled_end, scaled_end]
            ys = [
                annotation_base_y - ANNOTATION_HEIGHT,
                annotation_base_y + ANNOTATION_HEIGHT,
                annotation_base_y + ANNOTATION_HEIGHT,
                annotation_base_y - ANNOTATION_HEIGHT,
            ]

        # Update lowest annotation y
        if annotation_base_y - ANNOTATION_HEIGHT < lowest_annotation_y:
            lowest_annotation_y = annotation_base_y - ANNOTATION_HEIGHT

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

    lowest_annotation_y -= ANNOTATION_HEIGHT/2
    fig.line(
        [0, xmax],
        [lowest_annotation_y, lowest_annotation_y],
        line_width=2,
        color="black",
    )

    return fig, lowest_annotation_y

def plot_ssr(fig, ssr_df):
    height = 1
    ssr_shape = [[0, 0], [1, 0], [1, INNER_HEIGHT], [0, INNER_HEIGHT], [0, 0]]
    color = "green"

    ssr_source = {
        "x": [],
        "y": [],
        "color": [],
        "line_color": [],
        "name": [],
        "first_repeat": [],
        "mutation_rate": [],
        "line_width": [],
        "sequence": [],
    }

    outline_shape_left = [[0, 0], [1, 0],  [1, MARKER_HEIGHT], [0, MARKER_HEIGHT], [0, 0]]
    ssr_outline_source = {
        "x": [],
        "y": [],
        "color": [],
        "line_color": [],
        "name": [],
        "first_repeat": [],
        "mutation_rate": [],
        "line_width": [],
        "sequence": [],
        "count": [],
        "monomer": [],
    }

    for row in ssr_df.rows(named=True):
        left_pos = row["start"]
        repeat_len = row["repeat_len"]
        y_height = row["level"]*MARKER_HEIGHT*1.5
        glyph = [[point[0]*repeat_len+left_pos, point[1]+OUTLINE_PADDING_Y] for point in ssr_shape]

        ssr_source["x"].append([x for x, _ in glyph])
        ssr_source["y"].append([y+y_height+MARKER_HEIGHT*0.5 for _, y in glyph])
        ssr_source["color"].append("black")
        ssr_source["name"].append("ssr")
        ssr_source['first_repeat'].append(left_pos)
        ssr_source["mutation_rate"].append(row["mutation_rate"])
        ssr_source["line_width"].append(2)
        ssr_source["line_color"].append("black")
        ssr_source["sequence"].append(row["repeat"]*repeat_len)

        # Outline
        outline = [[point[0]*repeat_len+left_pos, point[1]+OUTLINE_PADDING_Y] for point in ssr_shape]

        outline_mods = [[-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y],
                        [OUTLINE_PADDING_X+repeat_len, -OUTLINE_PADDING_Y],
                        [OUTLINE_PADDING_X+repeat_len, OUTLINE_PADDING_Y],
                        [-OUTLINE_PADDING_X, OUTLINE_PADDING_Y],
                        [-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y]]
        for i in range(len(outline)):
            outline[i][0] = outline[i][0] + outline_mods[i][0]
            outline[i][1] = outline[i][1] + outline_mods[i][1] + y_height + MARKER_HEIGHT*0.5
        ssr_outline_source["x"].append([x for x, _ in outline])
        ssr_outline_source["y"].append([y for _, y in outline])
        ssr_outline_source["color"].append(None)
        ssr_outline_source["name"].append("SSR")
        ssr_outline_source['first_repeat'].append(left_pos)
        ssr_outline_source["mutation_rate"].append(row["mutation_rate"])
        ssr_outline_source["line_width"].append(3)
        ssr_outline_source["line_color"].append(color)
        ssr_outline_source["sequence"].append(row["repeat"]*row["count"])
        ssr_outline_source["monomer"].append(row["repeat"])
        ssr_outline_source['count'].append(row["count"])

    ssr_glyphs = fig.patches(
        "x",
        "y",
        color="color",
        source=ssr_source,
        alpha=0.5,
        line_color="line_color",
        line_width="line_width",
    )

    ssr_glyphs_outline = fig.patches(
        "x",
        "y",
        color="color",
        source=ssr_outline_source,
        alpha=0.5,
        line_color="line_color",
        line_width="line_width",
    )

    ssr_glyphs_hover = HoverTool(
        renderers=[ssr_glyphs_outline],
        tooltips=[
            ("Type", "@name"),
            ("Sequence", "@sequence"),
            ("Repeat Start", "@first_repeat"),
            ("Monomer", "@monomer"),
            ("Num. Repeats", "@count"),
            ("Mutation Rate", "@mutation_rate"),
        ],
    )

    fig.add_tools(ssr_glyphs_hover)

    return fig

def plot_srs(fig, ssr_df):
    height = 1
    rmd_shape = [[0, 0], [1, 0], [1, INNER_HEIGHT], [0, INNER_HEIGHT], [0, 0]]
    color = "orange"

    # Logic for exact position
    srs_source = {
        "x": [],
        "y": [],
        "color": [],
        "line_color": [],
        "name": [],
        "first_repeat": [],
        "second_repeat": [],
        "mutation_rate": [],
        "line_width": [],
        "sequence": [],
        "distance": [],
    }

    for row in ssr_df.rows(named=True):
        left_pos = row["first_repeat"]
        right_pos = row["second_repeat"]
        repeat_len = row["repeat_len"]
        y_height = row["level"]*MARKER_HEIGHT*1.5
        left_glyph = [[point[0]*repeat_len+left_pos, point[1]+OUTLINE_PADDING_Y] for point in rmd_shape]
        right_glyph = [[point[0]*repeat_len+right_pos, point[1]+OUTLINE_PADDING_Y] for point in rmd_shape]

        srs_source["x"].append([x for x, _ in left_glyph])
        srs_source["y"].append([y+y_height+MARKER_HEIGHT*0.5 for _, y in left_glyph])
        srs_source["x"].append([x for x, _ in right_glyph])
        srs_source["y"].append([y+y_height+ MARKER_HEIGHT*0.5 for _, y in right_glyph])
        for _ in range(2):
            srs_source["color"].append("black")
            srs_source["name"].append("srs")
            srs_source['first_repeat'].append(left_pos)
            srs_source['second_repeat'].append(right_pos)
            srs_source["mutation_rate"].append(row["mutation_rate"])
            srs_source["line_width"].append(2)
            srs_source["line_color"].append("black")
            srs_source["sequence"].append(row["repeat"])
            srs_source["distance"].append(row["distance"])

    srs_glyphs = fig.patches(
        "x",
        "y",
        color="color",
        source=srs_source,
        alpha=0.5,
        line_color="line_color",
        line_width="line_width",
    )


    # Logic for highlighting
    outline_shape_left = [[0, 0], [1, 0], [1, MARKER_HEIGHT/2], [1, MARKER_HEIGHT], [0, MARKER_HEIGHT], [0, 0]]
    outline_shape_right = [[0, 0], [1, 0], [1, MARKER_HEIGHT], [0, MARKER_HEIGHT], [0, MARKER_HEIGHT/2], [0, 0]]
    line_shape = [[1, MARKER_HEIGHT/2], [1,MARKER_HEIGHT/2]]
    srs_outline_source = {
        "x": [],
        "y": [],
        "color": [],
        "line_color": [],
        "name": [],
        "first_repeat": [],
        "second_repeat": [],
        "mutation_rate": [],
        "line_width": [],
        "sequence": [],
        "distance": [],
    }
    srs_line_source = {
        "x": [],
        "y": [],
        "color": [],
    }
    for row in ssr_df.rows(named=True):
        left_pos = row["first_repeat"]
        right_pos = row["second_repeat"]
        repeat_len = row["repeat_len"]
        y_height = row["level"]*MARKER_HEIGHT*1.5
        if row["distance"] > OUTLINE_PADDING_X*4: # Distant SRSs
            left_outline = [[point[0]+left_pos, point[1]+OUTLINE_PADDING_Y] for point in outline_shape_left]
            left_outline_mods = [[-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y],
                                 [OUTLINE_PADDING_X+repeat_len, -OUTLINE_PADDING_Y],
                                 [OUTLINE_PADDING_X*2+repeat_len, 0],
                                 [OUTLINE_PADDING_X+repeat_len, 0],
                                 [-OUTLINE_PADDING_X, 0],
                                 [-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y]]
            for i in range(len(left_outline)):
                left_outline[i][0] = left_outline[i][0] + left_outline_mods[i][0]
                left_outline[i][1] = left_outline[i][1] + left_outline_mods[i][1] + y_height + MARKER_HEIGHT*0.5

            right_outline = [[point[0]+right_pos, point[1]+OUTLINE_PADDING_Y] for point in outline_shape_right]
            right_outline_mods = [[-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y],
                                 [OUTLINE_PADDING_X+repeat_len, -OUTLINE_PADDING_Y],
                                 [OUTLINE_PADDING_X+repeat_len, 0],
                                 [-OUTLINE_PADDING_X+repeat_len, 0],
                                 [-OUTLINE_PADDING_X*2, 0],
                                 [-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y]]
            for i in range(len(right_outline)):
                right_outline[i][0] = right_outline[i][0] + right_outline_mods[i][0]
                right_outline[i][1] = right_outline[i][1] + right_outline_mods[i][1] + y_height + MARKER_HEIGHT*0.5

            srs_outline_source["x"].append([x for x, _ in left_outline])
            srs_outline_source["y"].append([y for _, y in left_outline])
            srs_outline_source["x"].append([x for x, _ in right_outline])
            srs_outline_source["y"].append([y for _, y in right_outline])
            for _ in range(2):
                srs_outline_source["color"].append(None)
                srs_outline_source["name"].append("srs")
                srs_outline_source['first_repeat'].append(left_pos)
                srs_outline_source['second_repeat'].append(right_pos)
                srs_outline_source["mutation_rate"].append(row["mutation_rate"])
                srs_outline_source["line_width"].append(3)
                srs_outline_source["line_color"].append(color)
                srs_outline_source["sequence"].append(row["repeat"])
                srs_outline_source["distance"].append(row["distance"])


            line = deepcopy(line_shape)
            line[0][0] = left_outline[2][0]
            line[1][0] = right_outline[4][0]
            srs_line_source["x"].append([x for x, _ in line])
            srs_line_source["y"].append([y+OUTLINE_PADDING_Y+ y_height + MARKER_HEIGHT*0.5 for _, y in line])
            srs_line_source["color"].append(color)

        else: # Near SRSs
            l = [[point[0]+left_pos, point[1]+OUTLINE_PADDING_Y] for point in rmd_shape]
            r = [[point[0]+right_pos, point[1]+OUTLINE_PADDING_Y] for point in rmd_shape]

            outline = [l[0], r[1], r[2], l[3], l[4]]
            outline_mods = [[-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y],
                            [OUTLINE_PADDING_X+repeat_len, -OUTLINE_PADDING_Y],
                            [OUTLINE_PADDING_X+repeat_len, OUTLINE_PADDING_Y],
                            [-OUTLINE_PADDING_X, OUTLINE_PADDING_Y],
                            [-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y]]
            for i in range(len(outline)):
                outline[i][0] = outline[i][0] + outline_mods[i][0]
                outline[i][1] = outline[i][1] + outline_mods[i][1] + y_height + MARKER_HEIGHT*0.5
            srs_outline_source["x"].append([x for x, _ in outline])
            srs_outline_source["y"].append([y for _, y in outline])
            srs_outline_source["color"].append(None)
            srs_outline_source["name"].append("SRS")
            srs_outline_source['first_repeat'].append(left_pos)
            srs_outline_source['second_repeat'].append(right_pos)
            srs_outline_source["mutation_rate"].append(row["mutation_rate"])
            srs_outline_source["line_width"].append(3)
            srs_outline_source["line_color"].append(color)
            srs_outline_source["sequence"].append(row["repeat"])
            srs_outline_source["distance"].append(row["distance"])

    srs_glyphs_outline = fig.patches(
        "x",
        "y",
        color="color",
        source=srs_outline_source,
        alpha=0.5,
        line_color="line_color",
        line_width="line_width",
    )

    srs_lines = fig.multi_line(
        xs="x",
        ys="y",
        source=srs_line_source,
        alpha=0.5,
        color="color",
        line_width=3,
    )

    srs_glyphs_hover = HoverTool(
        renderers=[srs_glyphs_outline],
        tooltips=[
            ("Type", "@name"),
            ("Sequence", "@sequence"),
            ("First Repeat (left)", "@first_repeat"),
            ("Second Repeat (right)", "@second_repeat"),
            ("Distance", "@distance"),
            ("Mutation Rate", "@mutation_rate"),
        ],
    )

    fig.add_tools(srs_glyphs_hover)

    return fig

def plot_rmd(fig, rmd_df):
    height = 1
    rmd_shape = [[0, 0], [1, 0], [1, INNER_HEIGHT], [0, INNER_HEIGHT], [0, 0]]
    color = "red"

    # Logic for exact position
    rmd_source = {
        "x": [],
        "y": [],
        "color": [],
        "line_color": [],
        "name": [],
        "first_repeat": [],
        "second_repeat": [],
        "mutation_rate": [],
        "line_width": [],
        "sequence": [],
        "distance": [],
    }

    for row in rmd_df.rows(named=True):
        left_pos = row["first_repeat"]
        right_pos = row["second_repeat"]
        repeat_len = row["repeat_len"]
        y_height = row["level"]*MARKER_HEIGHT*1.5
        left_glyph = [[point[0]*repeat_len+left_pos, point[1]+OUTLINE_PADDING_Y] for point in rmd_shape]
        right_glyph = [[point[0]*repeat_len+right_pos, point[1]+OUTLINE_PADDING_Y] for point in rmd_shape]

        rmd_source["x"].append([x for x, _ in left_glyph])
        rmd_source["y"].append([y+y_height+MARKER_HEIGHT*0.5 for _, y in left_glyph])
        rmd_source["x"].append([x for x, _ in right_glyph])
        rmd_source["y"].append([y+y_height+ MARKER_HEIGHT*0.5 for _, y in right_glyph])
        for _ in range(2):
            rmd_source["color"].append("black")
            rmd_source["name"].append("rmd")
            rmd_source['first_repeat'].append(left_pos)
            rmd_source['second_repeat'].append(right_pos)
            rmd_source["mutation_rate"].append(row["mutation_rate"])
            rmd_source["line_width"].append(2)
            rmd_source["line_color"].append("black")
            rmd_source["sequence"].append(row["repeat"])
            rmd_source["distance"].append(row["distance"])

    rmd_glyphs = fig.patches(
        "x",
        "y",
        color="color",
        source=rmd_source,
        alpha=0.5,
        line_color="line_color",
        line_width="line_width",
    )


    # Logic for highlighting
    outline_shape_left = [[0, 0], [1, 0], [1, MARKER_HEIGHT/2], [1, MARKER_HEIGHT], [0, MARKER_HEIGHT], [0, 0]]
    outline_shape_right = [[0, 0], [1, 0], [1, MARKER_HEIGHT], [0, MARKER_HEIGHT], [0, MARKER_HEIGHT/2], [0, 0]]
    line_shape = [[1, MARKER_HEIGHT/2], [1,MARKER_HEIGHT/2]]
    rmd_outline_source = {
        "x": [],
        "y": [],
        "color": [],
        "line_color": [],
        "name": [],
        "first_repeat": [],
        "second_repeat": [],
        "mutation_rate": [],
        "line_width": [],
        "sequence": [],
        "distance": [],
    }
    rmd_line_source = {
        "x": [],
        "y": [],
        "color": [],
    }
    for row in rmd_df.rows(named=True):
        left_pos = row["first_repeat"]
        right_pos = row["second_repeat"]
        repeat_len = row["repeat_len"]
        y_height = row["level"]*MARKER_HEIGHT*1.5
        if row["distance"] > OUTLINE_PADDING_X*4: # Distant rmds
            left_outline = [[point[0]+left_pos, point[1]+OUTLINE_PADDING_Y] for point in outline_shape_left]
            left_outline_mods = [[-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y],
                                 [OUTLINE_PADDING_X+repeat_len, -OUTLINE_PADDING_Y],
                                 [OUTLINE_PADDING_X*2+repeat_len, 0],
                                 [OUTLINE_PADDING_X+repeat_len, 0],
                                 [-OUTLINE_PADDING_X, 0],
                                 [-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y]]
            for i in range(len(left_outline)):
                left_outline[i][0] = left_outline[i][0] + left_outline_mods[i][0]
                left_outline[i][1] = left_outline[i][1] + left_outline_mods[i][1] + y_height + MARKER_HEIGHT*0.5

            right_outline = [[point[0]+right_pos, point[1]+OUTLINE_PADDING_Y] for point in outline_shape_right]
            right_outline_mods = [[-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y],
                                 [OUTLINE_PADDING_X+repeat_len, -OUTLINE_PADDING_Y],
                                 [OUTLINE_PADDING_X+repeat_len, 0],
                                 [-OUTLINE_PADDING_X, 0],
                                 [-OUTLINE_PADDING_X*2, 0],
                                 [-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y]]
            for i in range(len(right_outline)):
                right_outline[i][0] = right_outline[i][0] + right_outline_mods[i][0]
                right_outline[i][1] = right_outline[i][1] + right_outline_mods[i][1] + y_height + MARKER_HEIGHT*0.5

            rmd_outline_source["x"].append([x for x, _ in left_outline])
            rmd_outline_source["y"].append([y for _, y in left_outline])
            rmd_outline_source["x"].append([x for x, _ in right_outline])
            rmd_outline_source["y"].append([y for _, y in right_outline])
            for _ in range(2):
                rmd_outline_source["color"].append(None)
                rmd_outline_source["name"].append("rmd")
                rmd_outline_source['first_repeat'].append(left_pos)
                rmd_outline_source['second_repeat'].append(right_pos)
                rmd_outline_source["mutation_rate"].append(row["mutation_rate"])
                rmd_outline_source["line_width"].append(3)
                rmd_outline_source["line_color"].append(color)
                rmd_outline_source["sequence"].append(row["repeat"])
                rmd_outline_source["distance"].append(row["distance"])


            line = deepcopy(line_shape)
            line[0][0] = left_outline[2][0]
            line[1][0] = right_outline[4][0]
            rmd_line_source["x"].append([x for x, _ in line])
            rmd_line_source["y"].append([y+OUTLINE_PADDING_Y+ y_height + MARKER_HEIGHT*0.5 for _, y in line])
            rmd_line_source["color"].append(color)

        else: # Near rmds
            l = [[point[0]+left_pos, point[1]+OUTLINE_PADDING_Y] for point in rmd_shape]
            r = [[point[0]+right_pos, point[1]+OUTLINE_PADDING_Y] for point in rmd_shape]

            outline = [l[0], r[1], r[2], l[3], l[4]]
            outline_mods = [[-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y],
                            [OUTLINE_PADDING_X+repeat_len, -OUTLINE_PADDING_Y],
                            [OUTLINE_PADDING_X+repeat_len, OUTLINE_PADDING_Y],
                            [-OUTLINE_PADDING_X, OUTLINE_PADDING_Y],
                            [-OUTLINE_PADDING_X, -OUTLINE_PADDING_Y]]
            for i in range(len(outline)):
                outline[i][0] = outline[i][0] + outline_mods[i][0]
                outline[i][1] = outline[i][1] + outline_mods[i][1] + y_height + MARKER_HEIGHT*0.5
            rmd_outline_source["x"].append([x for x, _ in outline])
            rmd_outline_source["y"].append([y for _, y in outline])
            rmd_outline_source["color"].append(None)
            rmd_outline_source["name"].append("rmd")
            rmd_outline_source['first_repeat'].append(left_pos)
            rmd_outline_source['second_repeat'].append(right_pos)
            rmd_outline_source["mutation_rate"].append(row["mutation_rate"])
            rmd_outline_source["line_width"].append(3)
            rmd_outline_source["line_color"].append(color)
            rmd_outline_source["sequence"].append(row["repeat"])
            rmd_outline_source["distance"].append(row["distance"])

    rmd_glyphs_outline = fig.patches(
        "x",
        "y",
        color="color",
        source=rmd_outline_source,
        alpha=0.5,
        line_color="line_color",
        line_width="line_width",
    )

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
            ("First Repeat (left)", "@first_repeat"),
            ("Second Repeat (right)", "@second_repeat"),
            ("Distance", "@distance"),
            ("Mutation Rate", "@mutation_rate"),
        ],
    )

    fig.add_tools(rmd_glyphs_hover)

    return fig
