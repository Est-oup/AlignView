"""Presentation helpers shared by report builders."""
import xml.etree.ElementTree as ET


def figure_title(svg, fallback='Figure'):
    root = ET.fromstring(svg)
    title = root.find('.//{http://www.w3.org/2000/svg}text')
    return ''.join(title.itertext()) if title is not None else fallback
