#!/usr/bin/env python3
"""AlignView entry point and compatibility imports for existing Python callers."""
from av_cli import MODES, HELP_TEXT, build_parser, parse_args
from av_session import Session, show_file, text_segments, ascii_bar
from av_app import main, _apply_start, _select_start, show_tool_list

__all__ = ['main', 'Session', 'build_parser', 'parse_args', 'MODES', 'HELP_TEXT']

if __name__ == '__main__':
    raise SystemExit(main())
