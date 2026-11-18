"""Compatibility imports for reports created before AlignView 1.0."""
from av_figures import mapping_svg,density_svg


def build_report_v5(*args,**kwargs):
    from av_report import build_report
    return build_report(*args,**kwargs)
