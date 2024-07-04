from bokeh.models import ColumnDataSource, DataTable, TableColumn, CustomJS
from bokeh.models.widgets import Div
from bokeh.layouts import column
from bokeh.models.widgets.tables import CheckboxEditor


def generate_bokeh_table(df, name, callback=None) -> DataTable:
    # Generate a bokeh table from a polars dataframe
    column_names = df.columns
    table_name = name
    data = {column_name: [] for column_name in column_names}
    for datarow in df.iter_rows(named=True):
        for column_name in column_names:
            data[column_name].append(datarow[column_name])
    source = ColumnDataSource(data)
    columns = [
        TableColumn(field=column_name, title=column_name)
        for column_name in column_names
    ]
    name_div = Div(text=f"<h1>{name}</h1>")
    table = DataTable(
        source=source,
        columns=columns,
        name=table_name,
        width=1500,
        editable=True,
        selectable="checkbox",
    )
    table = column(name_div, table)
    if callback:
        callback = callback(source)
        source.selected.js_on_change("indices", callback)
    return table
