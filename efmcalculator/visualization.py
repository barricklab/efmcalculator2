import bokeh


def make_plot(seqrecord, ssr_dataframe=None, rmd_dataframe=None):
    view_format = "mirrored"
    # Set up plot

    annotation_scale_range = len(cryptresult.sequence)
    number_of_glyphs = math.floor(len(cryptresult.sequence) / 1000 * GLYPHS_PER_KB)

    # Set up the figure
    fig = figure(
        width=1500,
        height=750,
    )
    fig.xaxis.axis_label = "Position"
    fig.yaxis.axis_label = "Sequence of Concern"
    fig.yaxis.visible = False
    fig.ygrid.visible = False
    fig.xgrid.visible = False
    fig.toolbar.logo = None

    wigits = []
    tables = {}

    fig.add_layout(shownaxis, "left")

    fig.xaxis.axis_label_text_font_size = FONTSIZE
    fig.yaxis.axis_label_text_font_size = FONTSIZE
    shownaxis.axis_label_text_font_size = FONTSIZE
    fig.xaxis.axis_label_text_font = FONT
    fig.yaxis.axis_label_text_font = FONT
    shownaxis.axis_label_text_font = FONT
    fig.xaxis.axis_label_text_font_style = "bold"
    fig.yaxis.axis_label_text_font_style = "bold"
    shownaxis.axis_label_text_font_style = "bold"
    fig.xaxis.major_label_text_font_size = FONTSIZE
    fig.yaxis.major_label_text_font_size = FONTSIZE
    shownaxis.major_label_text_font_size = FONTSIZE
    fig.xaxis.major_label_text_font = FONT
    fig.yaxis.major_label_text_font = FONT
    shownaxis.major_label_text_font = FONT

    # set width to length of sequence
    fig.x_range = Range1d(start=0, end=len(cryptresult.sequence))

    # Draw a line for the top of the annotations
    fig.line(
        [0, annotation_scale_range],
        [0, 0],
        line_width=2,
        color="black",
        y_range_name="y_range2",
        x_range_name="x_range2",
    )

    if cryptresult.annotations:
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

        for genbank_annotation in cryptresult.annotations:
            arrow_depth = 100
            # @TODO: SCALE ANNOTATIONS

            scaled_start = (
                genbank_annotation.start
                * len(cryptresult.sequence)
                / annotation_scale_range
            )
            scaled_end = (
                genbank_annotation.end
                * len(cryptresult.sequence)
                / annotation_scale_range
            )

            if genbank_annotation.end - genbank_annotation.start < arrow_depth:
                arrow_depth = genbank_annotation.end - genbank_annotation.start
            annotation_base_y = (
                -1250.0 * genbank_annotation.nest_level + annotation_depth
            )

            if genbank_annotation.strand == 1:
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
            elif genbank_annotation.strand == -1:
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
            name = genbank_annotation.name

            # Update lowest annotation y
            if annotation_base_y - 500 < lowest_annotation_y:
                lowest_annotation_y = annotation_base_y - 500

            # Define the color
            color = genbank_annotation.color

            # Add the annotation to the dictionary

            genbank_dictionary["x"].append(xs)
            genbank_dictionary["y"].append(ys)
            genbank_dictionary["color"].append(color)
            genbank_dictionary["name"].append(name)
            genbank_dictionary["position"].append(
                f"{genbank_annotation.start}-{genbank_annotation.end}"
            )
            genbank_dictionary["strand"].append(genbank_annotation.strand)

        genbank_glyphs = fig.patches(
            "x",
            "y",
            color="color",
            source=genbank_dictionary,
            alpha=0.5,
            line_color="black",
            line_width=1,
            y_range_name="y_range2",
            x_range_name="x_range2",
        )
        genbank_glyphs_hover = HoverTool(
            renderers=[genbank_glyphs], tooltips=[("Name", "@name")]
        )
        fig.add_tools(genbank_glyphs_hover)

        # Draw a line below the annotations

        lowest_annotation_y -= 500
        fig.line(
            [0, annotation_scale_range],
            [lowest_annotation_y, lowest_annotation_y],
            line_width=2,
            color="black",
            y_range_name="y_range2",
            x_range_name="x_range2",
        )
        annotation_depth = lowest_annotation_y

    return fig


def generate_bokeh_table(datalist, name) -> DataTable:
    # Generate a bokeh table from a list of named tuples
    column_names = datalist[0]._fields
    table_name = datalist[0].__class__.__name__
    data = {column_name: [] for column_name in column_names}
    for datarow in datalist:
        for column_name in column_names:
            data[column_name].append(getattr(datarow, column_name))
    source = ColumnDataSource(data)
    columns = [
        TableColumn(field=column_name, title=column_name)
        for column_name in column_names
    ]
    name_div = Div(text=f"<h1>{name}</h1>")
    table = DataTable(
        source=source, columns=columns, name=table_name, width=1500, editable=True
    )
    table = column(name_div, table)
    return table
