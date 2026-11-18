"""Compatibility imports for reports created before AlignView 1.0."""
from av_figures import txt,shell,finish,hist_svg,bars_svg,scatter_svg,heatmap_svg,annotation_svg,plot_svg,INK,MUTED,BLUE,TEAL


def build_report_v4(*args,**kwargs):
    from av_report import build_report
    return build_report(*args,**kwargs)
