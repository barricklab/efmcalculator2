from bokeh.models import ColumnDataSource, DataTable, TableColumn
from bokeh.models.widgets import Div
from bokeh.layouts import column


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
