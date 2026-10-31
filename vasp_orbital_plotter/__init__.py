"""Orbital-projected band structures read directly from vasprun.xml.

Stage 1 provides the streaming parser (parser.py) and the "inspect" command.
Stage 2 adds selection and weights (projections.py), k-path distances and labels (kpath.py), plot defaults (options.py),
plots (plotting.py) and the plot command.
Stage 3 adds spin panels for collinear files, signed mx, my, mz plots, overlays, label renames, styles and large-file options.
The specification is in docs/orbital_plotter_spec.md.
"""

from .parser import (ProjectionsMissingError, VasprunData, VasprunError, VasprunLayoutError, VasprunParseError,
                     parse_vasprun)

__all__ = ["ProjectionsMissingError", "VasprunData", "VasprunError", "VasprunLayoutError", "VasprunParseError",
           "parse_vasprun"]
