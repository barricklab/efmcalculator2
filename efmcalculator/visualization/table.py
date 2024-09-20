from bokeh.models import ColumnDataSource, DataTable, TableColumn, CustomJS
from bokeh.models.widgets import Div
from bokeh.layouts import column
from bokeh.models.widgets.tables import CheckboxEditor, ScientificFormatter
import polars as pl


def generate_bokeh_table(df, name, callback=None) -> DataTable:
    # Generate a bokeh table from a polars dataframe
    column_names = df.columns
    table_name = name
    selected = df.filter(pl.col("show") == True).get_column("index").to_list()

    data = {column_name: [] for column_name in column_names}
    for datarow in df.iter_rows(named=True):
        for column_name in column_names:
            data[column_name].append(datarow[column_name])
    source = ColumnDataSource(data)
    columns = []
    for column_name in column_names:
        if column_name in ["index", "show"]:
            continue
        elif column_name == "mutation_rate":
            formatted_column = TableColumn(
                field=column_name,
                title=column_name,
                formatter=ScientificFormatter(precision=2),
            )
        else:
            formatted_column = TableColumn(field=column_name, title=column_name)
        columns.append(formatted_column)

    source.selected.indices = selected

    table = DataTable(
        source=source,
        columns=columns,
        name=table_name,
        width=500,
        editable=True,
        selectable="checkbox",
    )
    if callback:
        callback = callback(source)
        source.selected.js_on_change("indices", callback)
    return table

def generate_empty_table(name):
    table = Div(text = f"No {name} hotspots predicted", width=750)
    return table
