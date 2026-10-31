"""Orbital-projected band structures read directly from vasprun.xml.

Stage 1 provides the streaming parser (parser.py) and the "inspect" command.
Stage 2 adds selection and weights (projections.py), k-path distances and labels (kpath.py), plot defaults (options.py),
plots (plotting.py) and the plot command.
Spin panels, magnetisation plots and overlays are stage 3 of the specification in docs/orbital_plotter_spec.md.
"""

from .parser import (ProjectionsMissingError, VasprunData, VasprunError, VasprunLayoutError, VasprunParseError,
                     parse_vasprun)

__all__ = ["ProjectionsMissingError", "VasprunData", "VasprunError", "VasprunLayoutError", "VasprunParseError",
           "parse_vasprun"]
