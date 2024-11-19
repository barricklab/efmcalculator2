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
                title=label_rename(column_name),
                formatter=ScientificFormatter(precision=2),
            )
        else:
            formatted_column = TableColumn(field=column_name, title=label_rename(column_name))
        columns.append(formatted_column)

    source.selected.indices = selected

    table = DataTable(
        source=source,
        columns=columns,
        name=table_name,
        width=500,
        editable=True,
        selectable="checkbox",
        index_position = None
    )
    if callback:
        callback = callback(source)
        source.selected.js_on_change("indices", callback)
    return table

def generate_empty_table(name):
    table = Div(text = f"No {name} hotspots predicted", width=750)
    return table

def label_rename(title):
    result = title.replace("_", " ")
        
    if len(title) > 0:
        result = result[0].upper() + result[1: len(result)]
    
    if result.find(" ") != -1 and result.find(" ") != (len(result) - 1):
        space_index = result.find(" ")
        result = result[0 : space_index + 1] + result[space_index + 1].upper() + result[space_index + 2 : len(result)]
        
    return result

def generate_nerfed_bokeh_table(polarTable):
    if isinstance(polarTable, pl.DataFrame):
        polarTable = polarTable.to_pandas()
        polarTable["mutation_rate"] = polarTable["mutation_rate"].apply(lambda x: f"{x:.2e}")
    
    source = ColumnDataSource(polarTable)
    columns = [TableColumn(field=col, title=label_rename(col)) for col in polarTable.columns]
    data_table = DataTable(source=source, columns=columns, width=500, index_position = None)
    
    return data_table