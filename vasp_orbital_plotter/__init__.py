"""Orbital-projected band structures read directly from vasprun.xml.

Stage 1 (this version) provides the streaming parser in vasp_orbital_plotter.parser and the "inspect" command.
Projections, k-path reconstruction and plotting are stages 2 and 3 of the specification in docs/orbital_plotter_spec.md.
"""

from .parser import (ProjectionsMissingError, VasprunData, VasprunError, VasprunLayoutError, VasprunParseError,
                     parse_vasprun)

__all__ = ["ProjectionsMissingError", "VasprunData", "VasprunError", "VasprunLayoutError", "VasprunParseError",
           "parse_vasprun"]
