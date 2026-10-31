"""Fat-band and colour-projected band plots with Matplotlib (stages 2 and 3 of the specification).

The figure is built on a matplotlib.figure.Figure with the Agg canvas, so no window and no global backend change is needed.
Bands are drawn per path branch, so separate branches are never joined.
Fat-band mode draws every band as a thin grey line and puts one marker per state, with area marker_scale times the weight.
Colour mode draws each band in colour; every half-step next to a computed state takes that state's colour, so no colours are blended.
A figure has one or two panels (spin up and down) sharing the energy axis, each with one or more series (overlays).
Signed weights (mx, my, mz) use a diverging map with a grey middle centred on zero, or red and blue markers sized by |value|.
With several series, colour mode gives each state the colour of its dominant series, faded to grey as that weight falls to zero.
"""

from __future__ import annotations

import textwrap
from dataclasses import dataclass, replace
from functools import partial
from itertools import repeat
from pathlib import Path
from typing import Sequence

import matplotlib
import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Colormap, LinearSegmentedColormap, ListedColormap, Normalize, TwoSlopeNorm, to_rgba
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from altermagagent.bands import mpl_label, segments

from .kpath import KPath
from .options import DEFAULT_WINDOW, PANEL_SIZE, SERIES_COLOURS, OptionsError, PlotOptions
from .parser import VasprunData

BAND_GREY = "#9a9a9a"
# Marker colour of a single selection: the first series colour (SERIES_COLOURS lives in options.py, free of Matplotlib).
MARKER_COLOUR = SERIES_COLOURS[0]
# Signed fat bands: positive and negative values, the two poles of RdBu.
POSITIVE_COLOUR, NEGATIVE_COLOUR = "#b2182b", "#2166ac"
# Neutral middle of the signed colour map and the faded end of dominant-series colours: grey, so zero stays visible on white.
NEUTRAL_GREY = "#c4c4c4"
CAPTION_GREY = "#555555"
# Fraction cut from the light end of the colour map, so states of zero weight stay visible on white paper.
CMAP_LIGHT_CUT = 0.2


@dataclass
class Series:
    """One projection drawn on a panel: weights of shape (kpoint, band), its legend label and its colour."""

    weights: np.ndarray
    label: str
    colour: str = MARKER_COLOUR


@dataclass
class Panel:
    """One band panel: energies of shape (kpoint, band), its series, and a heading such as "spin up" (None for none)."""

    energies: np.ndarray
    series: list[Series]
    heading: str | None = None


# Fermi energy and the energy window in plotted units: the window centre is E_F for absolute plots and 0 for E - E_F.
def resolve_window(data: VasprunData, emin: float | None, emax: float | None, absolute: bool) -> tuple[float, float, float]:
    if data.efermi is None:
        raise OptionsError("the file has no Fermi energy; give one with --efermi")
    centre = data.efermi if absolute else 0.0
    return (data.efermi, centre - DEFAULT_WINDOW if emin is None else emin,
            centre + DEFAULT_WINDOW if emax is None else emax)


# Colour map with its lightest CMAP_LIGHT_CUT removed, so the lowest weights are pale but not white.
def visible_cmap(name: str) -> ListedColormap:
    base = matplotlib.colormaps[name]
    return ListedColormap(base(np.linspace(CMAP_LIGHT_CUT, 1.0, 256)), name=f"{name}_visible")


# Diverging colour map whose middle is neutral grey, so states with zero magnetisation stay visible on white paper.
def grey_centred_cmap(name: str) -> Colormap:
    base = matplotlib.colormaps[name]
    return LinearSegmentedColormap.from_list(f"{name}_grey_centre", [base(0.0), base(0.2), NEUTRAL_GREY, base(0.8), base(1.0)])


# Points of every band on one branch, shape (kpoint, band, 2), from distances (kpoint,) and energies (kpoint, band).
def _points(xs: np.ndarray, es: np.ndarray) -> np.ndarray:
    return np.stack([np.broadcast_to(xs[:, None], es.shape), es], axis=-1)


# Line pieces of every band on every branch: whole steps, or with values half steps split at the midpoints.
def _branch_pieces(kpath: KPath, es: np.ndarray, values: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray | None]:
    """With values (kpoint, band, ...), each half step also returns the value of the state it touches; otherwise that is None."""
    pieces, touched = [], []
    value_parts = segments(kpath.x, values, kpath.branches) if values is not None else repeat((None, None))
    for (xs, ys), (_, v) in zip(segments(kpath.x, es, kpath.branches), value_parts):
        p = _points(xs, ys)
        if v is None:
            pieces.append(np.stack([p[:-1], p[1:]], axis=2).reshape(-1, 2, 2))
            continue
        mid = (p[:-1] + p[1:]) / 2
        pieces += [np.stack([p[:-1], mid], axis=2).reshape(-1, 2, 2), np.stack([mid, p[1:]], axis=2).reshape(-1, 2, 2)]
        touched += [v[:-1].reshape(-1, *v.shape[2:]), v[1:].reshape(-1, *v.shape[2:])]
    return np.concatenate(pieces), np.concatenate(touched) if touched else None


# Thin grey line of every band on every branch, drawn under the fat-band markers.
def _grey_bands(ax: Axes, kpath: KPath, es: np.ndarray, options: PlotOptions) -> None:
    pieces, _ = _branch_pieces(kpath, es)
    ax.add_collection(LineCollection(pieces, colors=BAND_GREY, linewidths=options.linewidth * 0.6, zorder=1,
                                     rasterized=options.rasterize))


# Fat-band mode: grey band lines, then markers of the states inside the window; signed series get a colour per sign.
def _draw_fatband(ax: Axes, kpath: KPath, panel: Panel, options: PlotOptions, signed: bool) -> None:
    es = panel.energies
    _grey_bands(ax, kpath, es, options)
    inside = (es >= options.emin) & (es <= options.emax)
    xs = np.broadcast_to(kpath.x[:, None], es.shape)
    for series in panel.series:
        ws = series.weights
        parts = ((ws > 0, POSITIVE_COLOUR), (ws < 0, NEGATIVE_COLOUR)) if signed else ((ws > 0, series.colour),)
        for keep, colour in parts:
            show = inside & keep
            ax.scatter(xs[show], es[show], s=options.marker_scale * np.abs(ws[show]), color=colour, alpha=options.alpha,
                       linewidths=0, zorder=2, rasterized=options.rasterize)


# Colour-mode norm shared by all panels: 0 to the window maximum, or symmetric about zero for signed weights.
def _color_norm(panels: Sequence[Panel], options: PlotOptions, signed: bool) -> Normalize:
    visible = [s.weights[(p.energies >= options.emin) & (p.energies <= options.emax)] for p in panels for s in p.series]
    values = np.concatenate(visible) if visible else np.zeros(0)
    if signed:
        top = float(np.abs(values).max()) if values.size else 0.0
        top = top if top > 0 else 1.0
        low = -top if options.cmin is None else options.cmin
        high = top if options.cmax is None else options.cmax
        if not low < 0 < high:
            raise OptionsError(f"signed colour limits must lie on both sides of zero, not {low:g} to {high:g}")
        return TwoSlopeNorm(vcenter=0.0, vmin=low, vmax=high)
    top = float(values.max()) if values.size else 1.0
    return Normalize(vmin=0.0 if options.cmin is None else options.cmin,
                     vmax=(top if top > 0 else 1.0) if options.cmax is None else options.cmax)


# Colour mode with one series: every band as half-steps coloured by the weight of the state they touch.
def _draw_color(ax: Axes, kpath: KPath, panel: Panel, options: PlotOptions, norm: Normalize, cmap: Colormap) -> None:
    pieces, values = _branch_pieces(kpath, panel.energies, panel.series[0].weights)
    lines = LineCollection(pieces, cmap=cmap, norm=norm, linewidths=options.linewidth * 2.2, zorder=2, capstyle="round",
                           rasterized=options.rasterize)
    lines.set_array(values)
    ax.add_collection(lines)


# Colour mode with several series: each state takes its dominant series' colour, faded to grey as that weight falls to zero.
def _draw_dominant(ax: Axes, kpath: KPath, panel: Panel, options: PlotOptions, norm: Normalize) -> None:
    stack = np.stack([np.clip(s.weights, 0.0, None) for s in panel.series])
    dominant, strength = stack.argmax(axis=0), stack.max(axis=0)
    fade = np.clip(norm(strength), 0.0, 1.0)[..., None]
    colours = np.array([to_rgba(s.colour)[:3] for s in panel.series])[dominant]
    rgb = np.array(to_rgba(NEUTRAL_GREY)[:3]) * (1 - fade) + colours * fade
    pieces, touched = _branch_pieces(kpath, panel.energies, rgb)
    ax.add_collection(LineCollection(pieces, colors=touched, linewidths=options.linewidth * 2.2, zorder=2, capstyle="round",
                                     rasterized=options.rasterize))


# Legend handles: one marker or line per series, one marker per sign for signed fat bands, none for signed colour mode.
def _legend_handles(panel: Panel, options: PlotOptions, signed: bool, quantity: str) -> list[Line2D]:
    """Signed colour mode has no entry: a flat-colour line would hide the sign, and the colour bar already names the quantity."""
    size = options.font_size * 0.6
    if signed and options.mode == "color":
        return []
    if signed and options.mode == "fatband":
        return [Line2D([], [], linestyle="", marker="o", markersize=size, color=colour, alpha=options.alpha,
                       label=f"{panel.series[0].label}: {quantity} {sign} 0")
                for colour, sign in ((POSITIVE_COLOUR, ">"), (NEGATIVE_COLOUR, "<"))]
    if options.mode == "fatband":
        return [Line2D([], [], linestyle="", marker="o", markersize=size, color=s.colour, alpha=options.alpha, label=s.label)
                for s in panel.series]
    return [Line2D([], [], color=s.colour, linewidth=options.linewidth * 2.2, label=s.label) for s in panel.series]


# Draw and save a figure of one or more panels sharing the energy axis; when no band reaches the window, add to notes.
def plot_panels(panels: Sequence[Panel], kpath: KPath, output: Path, options: PlotOptions, zero: float = 0.0,
                notes: list[str] | None = None, signed: bool = False, quantity: str = "projection weight",
                caption: str | None = None) -> Path:
    """Energies are in plotted units (E - E_F unless --absolute) and zero is the Fermi level in those units.

    Panels are drawn as given: the caller cuts them to the bands that reach the window, so an empty panel means none does.
    Without options.figsize the figure is PANEL_SIZE wide per panel; this is the only place that default is computed.
    """
    output = Path(output)
    options.validate(output)
    if options.emin is None or options.emax is None:
        raise ValueError("the energy window is not set; resolve it with resolve_window first")
    if not panels:
        raise ValueError("no panels to draw")
    for i, panel in enumerate(panels, 1):
        if not panel.series:
            raise ValueError(f"panel {i} has no series; give it at least one Series of weights")
        if panel.energies.shape[0] != kpath.x.shape[0]:
            raise ValueError(f"panel {i} has energies for {panel.energies.shape[0]} k-points, "
                             f"but the k-path has {kpath.x.shape[0]}")
        for series in panel.series:
            if panel.energies.shape != series.weights.shape:
                raise ValueError(f"panel {i}: energies {panel.energies.shape} and weights {series.weights.shape} "
                                 f"of {series.label!r} do not match")
    if options.figsize is None:
        options = replace(options, figsize=(PANEL_SIZE[0] * len(panels), PANEL_SIZE[1]))
    options = options.styled()
    with matplotlib.rc_context(options.rc()):
        fig = _figure(panels, kpath, options, zero, notes, signed, quantity, caption)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=options.dpi)
    return output


# Build the figure: panels, colour bar, legend, title and caption; options.figsize is the final size set by plot_panels.
def _figure(panels: Sequence[Panel], kpath: KPath, options: PlotOptions, zero: float, notes: list[str] | None,
            signed: bool, quantity: str, caption: str | None) -> Figure:
    if options.figsize is None:
        raise ValueError("the figure size is not set; plot_panels sets its default")
    figsize = options.figsize
    fig = Figure(figsize=figsize, layout="constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, len(panels), sharey=True, squeeze=False)[0]
    drawn = [p.energies.shape[1] > 0 for p in panels]
    several = any(len(p.series) > 1 for p in panels)
    legend_wanted = (options.legend or several or (signed and options.mode == "fatband")) and any(drawn)
    mappable = None
    if options.mode == "fatband":
        draw = partial(_draw_fatband, kpath=kpath, options=options, signed=signed)
    else:
        norm = _color_norm([p for p, d in zip(panels, drawn) if d], options, signed) if any(drawn) else None
        if several:
            draw = partial(_draw_dominant, kpath=kpath, options=options, norm=norm)
        else:
            cmap = grey_centred_cmap(options.cmap_signed) if signed else visible_cmap(options.cmap)
            draw = partial(_draw_color, kpath=kpath, options=options, norm=norm, cmap=cmap)
            mappable = ScalarMappable(norm=norm, cmap=cmap) if norm is not None else None
    for ax, panel, ok in zip(axes, panels, drawn):
        if ok:
            draw(ax, panel=panel)
    if not any(drawn) and notes is not None:
        unit = "eV" if options.absolute else "eV relative to the Fermi level"
        notes.append(f"no bands in the window {options.emin:g} to {options.emax:g} {unit}; "
                     "the figure has empty axes and no legend")
    for i, (ax, panel) in enumerate(zip(axes, panels)):
        _decorate(ax, kpath, options, zero, first=i == 0, heading=panel.heading if len(panels) > 1 else None)
    if options.title and len(panels) > 1:
        fig.suptitle(options.title, fontsize=options.font_size)
    elif options.title:
        axes[0].set_title(options.title, fontsize=options.font_size)
    if mappable is not None:
        bar = fig.colorbar(mappable, ax=list(axes), pad=0.02)
        what = f"{quantity} projection (sign convention unverified)" if signed else (
            "normalised weight" if options.normalize else "projection weight")
        bar.set_label(what, fontsize=options.font_size)
        bar.ax.tick_params(labelsize=options.font_size * 0.85, direction="in")
    handles = _legend_handles(panels[0], options, signed, quantity) if legend_wanted else []
    if handles:
        axes[0].legend(handles=handles, loc="best", fontsize=options.font_size * 0.8, frameon=True, framealpha=0.85,
                       edgecolor="none")
    if caption:
        size = options.font_size * 0.7
        # Wrapped to the figure width, taking an average character as 0.55 of the font size, so no caption is cut off.
        width = max(20, int(figsize[0] * 72 / (0.55 * size)))
        fig.supxlabel(textwrap.fill(caption, width), fontsize=size, color=CAPTION_GREY)
    return fig


# Axes limits, high-symmetry lines and labels, zero line at the Fermi level, panel heading and the frame style.
def _decorate(ax: Axes, kpath: KPath, options: PlotOptions, zero: float, first: bool = True,
              heading: str | None = None) -> None:
    x = kpath.x
    ax.set_xlim(float(x[0]), float(x[-1]))
    ax.set_ylim(options.emin, options.emax)
    positions = sorted({float(x[s]) for s in kpath.segment_starts} | {p for p, _ in kpath.ticks})
    for position in positions:
        if x[0] < position < x[-1]:
            ax.axvline(position, color="#555555", linewidth=0.6, zorder=0)
    labelled = [(p, lab) for p, lab in kpath.ticks if lab]
    ax.set_xticks([p for p, _ in labelled], [mpl_label(lab) for _, lab in labelled])
    if options.emin < zero < options.emax:
        ax.axhline(zero, color="#333333", linewidth=0.6, linestyle="--", zorder=0)
    if first:
        ax.set_ylabel("Energy (eV)" if options.absolute else r"$E - E_\mathrm{F}$ (eV)", fontsize=options.font_size)
    ax.tick_params(labelsize=options.font_size * 0.9, direction="in", top=False, right=not options.despine,
                   labelleft=first)
    ax.tick_params(axis="x", length=0)
    if options.despine:
        ax.spines[["top", "right"]].set_visible(False)
    if heading:
        ax.set_title(heading, fontsize=options.font_size * 0.9)
