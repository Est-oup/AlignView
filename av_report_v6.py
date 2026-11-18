"""Compatibility imports for reports created before AlignView 1.0."""
from av_figures import tracks_svg,frame_tracks_svg,atlas_additions


def build_report_v6(*args,**kwargs):
    from av_report import build_report
    return build_report(*args,**kwargs)
