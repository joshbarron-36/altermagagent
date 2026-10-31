"""Plot modes, output formats, the PlotOptions defaults and the light shared helpers, kept free of Matplotlib.

The command line reads its defaults from PlotOptions, so each default is written once, here.
OptionsError and KPathError live here, so the command line can catch them without importing the plotting or k-path modules.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

MODES = ("fatband", "color")
FORMATS = (".png", ".svg", ".pdf")
# Half width in eV of the default energy window around the Fermi level (or around 0 for E - E_F plots).
DEFAULT_WINDOW = 3.0


class OptionsError(ValueError):
    """A plot option has a value the figure cannot be drawn with."""


class KPathError(ValueError):
    """A label option or k-path input is written wrongly or does not fit the k-point list."""


# Split a comma list such as "1,3-4,Mn" or "dxy, eg" into stripped, non-empty items.
def split_list(text: str | None) -> list[str]:
    return [item.strip() for item in (text or "").split(",") if item.strip()]


@dataclass
class PlotOptions:
    """Look of one figure; every field has a publication default, and None for emin or emax means the default window."""

    mode: str = "fatband"
    emin: float | None = None
    emax: float | None = None
    absolute: bool = False
    normalize: bool = False
    figsize: tuple[float, float] = (6.0, 4.5)
    dpi: int = 300
    font_size: float = 12.0
    linewidth: float = 0.8
    marker_scale: float = 60.0
    alpha: float = 0.6
    cmap: str = "Blues"
    cmin: float | None = None
    cmax: float | None = None
    title: str | None = None
    legend: bool = False

    # Raise OptionsError for values Matplotlib would reject with a traceback; output, when given, must have a known suffix.
    def validate(self, output: Path | None = None) -> None:
        if self.mode not in MODES:
            raise OptionsError(f"mode {self.mode!r} is not one of {', '.join(MODES)}")
        if output is not None and Path(output).suffix.lower() not in FORMATS:
            raise OptionsError(f"output {Path(output).name!r} must end in {', '.join(FORMATS)}")
        if self.emin is not None and self.emax is not None and not self.emin < self.emax:
            raise OptionsError(f"the lower energy limit {self.emin:g} must be below the upper limit {self.emax:g}")
        if self.figsize[0] <= 0 or self.figsize[1] <= 0 or self.dpi <= 0:
            raise OptionsError("figure size and dpi must be positive")
        if self.marker_scale < 0 or not 0 <= self.alpha <= 1:
            raise OptionsError("the marker scale must be 0 or more and the opacity between 0 and 1")
        if self.font_size <= 0:
            raise OptionsError(f"the font size must be positive, not {self.font_size:g}")
        if self.linewidth < 0:
            raise OptionsError(f"the line width must be 0 or more, not {self.linewidth:g}")
        if self.cmin is not None and self.cmax is not None and not self.cmin < self.cmax:
            raise OptionsError(f"the lower colour limit {self.cmin:g} must be below the upper limit {self.cmax:g}")
        import matplotlib  # only here, so building the command line never loads Matplotlib
        if self.cmap not in matplotlib.colormaps:
            raise OptionsError(f"colour map {self.cmap!r} is not known to Matplotlib")
