"""Plot modes, output formats, the PlotOptions defaults and the light shared helpers, kept free of Matplotlib.

The command line reads its defaults from PlotOptions, so each default is written once, here.
Styles scale the base sizes, so an explicit --font-size keeps its meaning under --style talk.
OptionsError and KPathError live here, so the command line can catch them without importing the plotting or k-path modules.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

MODES = ("fatband", "color")
FORMATS = (".png", ".svg", ".pdf")
# Half width in eV of the default energy window around the Fermi level (or around 0 for E - E_F plots).
DEFAULT_WINDOW = 3.0
# Width and height in inches of one panel; a figure of several panels side by side is that many panels wide.
PANEL_SIZE = (6.0, 4.5)
# rcParams of every style: inward ticks everywhere, and TrueType fonts in PDF and PS so text stays editable.
_BASE_RC = {"xtick.direction": "in", "ytick.direction": "in", "pdf.fonttype": 42, "ps.fonttype": 42}
# Each --style: factors on (font size, line width, marker area) and its rcParams; paper is the stage 2 look.
STYLES: dict[str, dict[str, Any]] = {"paper": {"scale": (1.0, 1.0, 1.0), "rc": _BASE_RC},
                                     "talk": {"scale": (1.4, 1.5, 2.0), "rc": _BASE_RC}}
# Series colours of an overlay, in this fixed order (stage 2 blue first); checked for colour-vision-deficiency separation.
SERIES_COLOURS = ("#1f4e9c", "#d55e00", "#009e73", "#e69f00", "#cc79a7", "#56b4e9")
# Most series one overlay may have: one per colour.
MAX_SERIES = len(SERIES_COLOURS)


class OptionsError(ValueError):
    """A plot option has a value the figure cannot be drawn with."""


class KPathError(ValueError):
    """A label option or k-path input is written wrongly or does not fit the k-point list."""


# Split a comma list such as "1,3-4,Mn" or "dxy, eg" into stripped, non-empty items.
def split_list(text: str | None) -> list[str]:
    return [item.strip() for item in (text or "").split(",") if item.strip()]


@dataclass
class PlotOptions:
    """Look of one figure; every field has a publication default, and None for emin or emax means the default window.

    figsize None means PANEL_SIZE per panel. rcparams holds the parsed --rcparams-file once validate has read it; it has no flag.
    """

    mode: str = "fatband"
    emin: float | None = None
    emax: float | None = None
    absolute: bool = False
    normalize: bool = False
    figsize: tuple[float, float] | None = None
    dpi: int = 300
    font_size: float = 12.0
    linewidth: float = 0.8
    marker_scale: float = 60.0
    alpha: float = 0.6
    cmap: str = "Blues"
    cmin: float | None = None
    cmax: float | None = None
    cmap_signed: str = "RdBu_r"
    title: str | None = None
    legend: bool = False
    style: str = "paper"
    rcparams_file: str | None = None
    despine: bool = False
    rasterize: bool = False
    rcparams: dict[str, Any] | None = field(default=None, compare=False, repr=False, metadata={"flag": False})

    # Copy with font size, line width and marker area multiplied by the factors of the chosen style.
    def styled(self) -> "PlotOptions":
        font, line, marker = STYLES[self.style]["scale"]
        return replace(self, font_size=self.font_size * font, linewidth=self.linewidth * line,
                       marker_scale=self.marker_scale * marker)

    # Raise OptionsError for values Matplotlib would reject; output needs a known suffix; the rcParams file is read once here.
    def validate(self, output: Path | None = None) -> None:
        if self.mode not in MODES:
            raise OptionsError(f"mode {self.mode!r} is not one of {', '.join(MODES)}")
        if output is not None and Path(output).suffix.lower() not in FORMATS:
            raise OptionsError(f"output {Path(output).name!r} must end in {', '.join(FORMATS)}")
        if self.emin is not None and self.emax is not None and not self.emin < self.emax:
            raise OptionsError(f"the lower energy limit {self.emin:g} must be below the upper limit {self.emax:g}")
        if (self.figsize is not None and min(self.figsize) <= 0) or self.dpi <= 0:
            raise OptionsError("figure size and dpi must be positive")
        if self.marker_scale < 0 or not 0 <= self.alpha <= 1:
            raise OptionsError("the marker scale must be 0 or more and the opacity between 0 and 1")
        if self.font_size <= 0:
            raise OptionsError(f"the font size must be positive, not {self.font_size:g}")
        if self.linewidth < 0:
            raise OptionsError(f"the line width must be 0 or more, not {self.linewidth:g}")
        if self.cmin is not None and self.cmax is not None and not self.cmin < self.cmax:
            raise OptionsError(f"the lower colour limit {self.cmin:g} must be below the upper limit {self.cmax:g}")
        if self.style not in STYLES:
            raise OptionsError(f"style {self.style!r} is not one of {', '.join(STYLES)}")
        import matplotlib  # only here, so building the command line never loads Matplotlib
        for name in (self.cmap, self.cmap_signed):
            if name not in matplotlib.colormaps:
                raise OptionsError(f"colour map {name!r} is not known to Matplotlib")
        if self.rcparams_file is not None and self.rcparams is None:
            self.rcparams = load_rcparams(self.rcparams_file)

    # rcParams of the figure: those of the style, then those of the rcParams file read by validate.
    def rc(self) -> dict[str, Any]:
        return {**STYLES[self.style]["rc"], **(self.rcparams or {})}


# Read a JSON object of Matplotlib rcParams and check every key and value; problems raise OptionsError.
def load_rcparams(path: str | Path) -> dict[str, Any]:
    path = Path(path).expanduser()
    try:
        values = json.loads(path.read_text())
    except OSError as error:
        raise OptionsError(f"rcParams file {path} cannot be read ({error.strerror or error})") from None
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise OptionsError(f"rcParams file {path} is not valid JSON ({error})") from None
    if not isinstance(values, dict):
        raise OptionsError(f"rcParams file {path} must hold one JSON object, such as {{\"font.family\": \"serif\"}}")
    import matplotlib
    try:
        matplotlib.RcParams(values)
    except (KeyError, ValueError) as error:
        message = error.args[0] if error.args else error
        raise OptionsError(f"rcParams file {path}: {message}") from None
    return values
