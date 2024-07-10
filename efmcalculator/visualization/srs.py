from .rmd import _draw_rmd_logic


def draw_srs(fig, rmd_df):
    color = "orange"
    return _draw_rmd_logic(fig, rmd_df, color, "SRS")
