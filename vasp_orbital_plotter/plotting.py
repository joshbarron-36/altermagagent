"""Fat-band and colour-projected band plots with Matplotlib (stage 2 of the specification).

The figure is built on a matplotlib.figure.Figure with the Agg canvas, so no window and no global backend change is needed.
Bands are drawn per path branch, so separate branches are never joined.
Fat-band mode draws every band as a thin grey line and puts one marker per state, with area marker_scale times the weight.
Colour mode draws each band in colour; every half-step next to a computed state takes that state's colour, so no colours are blended.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import matplotlib
import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.figure import Figure

from altermagagent.bands import mpl_label, segments

from .kpath import KPath
from .options import DEFAULT_WINDOW, OptionsError, PlotOptions
from .parser import VasprunData

BAND_GREY = "#9a9a9a"
MARKER_COLOUR = "#1f4e9c"
# Fraction cut from the light end of the colour map, so states of zero weight stay visible on white paper.
CMAP_LIGHT_CUT = 0.2


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


# Points of every band on one branch, shape (kpoint, band, 2), from distances (kpoint,) and energies (kpoint, band).
def _points(xs: np.ndarray, es: np.ndarray) -> np.ndarray:
    return np.stack([np.broadcast_to(xs[:, None], es.shape), es], axis=-1)


# Fat-band mode: grey band lines, then one scatter of the states inside the window with a positive weight.
def _draw_fatband(ax: Axes, kpath: KPath, es: np.ndarray, ws: np.ndarray, options: PlotOptions, label: str) -> None:
    grey = [np.stack([p[:-1], p[1:]], axis=2).reshape(-1, 2, 2)
            for p in (_points(xs, ys) for xs, ys in segments(kpath.x, es, kpath.branches))]
    ax.add_collection(LineCollection(np.concatenate(grey), colors=BAND_GREY, linewidths=options.linewidth * 0.6, zorder=1))
    show = (es >= options.emin) & (es <= options.emax) & (ws > 0)
    xs = np.broadcast_to(kpath.x[:, None], es.shape)
    ax.scatter(xs[show], es[show], s=options.marker_scale * ws[show], color=MARKER_COLOUR, alpha=options.alpha,
               linewidths=0, zorder=2, label=label)


# Colour mode: every band as half-steps coloured by the weight of the state they touch, with a colour bar.
def _draw_color(ax: Axes, fig: Figure, kpath: KPath, es: np.ndarray, ws: np.ndarray, options: PlotOptions,
                label: str) -> None:
    pieces, values = [], []
    for (xs, ys), (_, weights) in zip(segments(kpath.x, es, kpath.branches), segments(kpath.x, ws, kpath.branches)):
        p = _points(xs, ys)
        mid = (p[:-1] + p[1:]) / 2
        pieces += [np.stack([p[:-1], mid], axis=2).reshape(-1, 2, 2), np.stack([mid, p[1:]], axis=2).reshape(-1, 2, 2)]
        values += [weights[:-1].ravel(), weights[1:].ravel()]
    visible = (es >= options.emin) & (es <= options.emax)
    top = float(ws[visible].max()) if visible.any() else 1.0
    norm = Normalize(vmin=0.0 if options.cmin is None else options.cmin,
                     vmax=(top if top > 0 else 1.0) if options.cmax is None else options.cmax)
    lines = LineCollection(np.concatenate(pieces), cmap=visible_cmap(options.cmap), norm=norm,
                           linewidths=options.linewidth * 2.2, zorder=2, capstyle="round", label=label)
    lines.set_array(np.concatenate(values))
    ax.add_collection(lines)
    bar = fig.colorbar(lines, ax=ax, pad=0.02)
    bar.set_label("normalised weight" if options.normalize else "projection weight", fontsize=options.font_size)
    bar.ax.tick_params(labelsize=options.font_size * 0.85)


# Drawer of each mode, called as draw(ax, fig, kpath, energies, weights, options, label); fat bands need no figure.
DRAW: dict[str, Callable[..., None]] = {
    "fatband": lambda ax, fig, *rest: _draw_fatband(ax, *rest),
    "color": _draw_color,
}


# Draw and save the figure from (kpoint, band) arrays; when no band reaches the window, save empty axes and add to notes.
def plot_bands(energies: np.ndarray, weights: np.ndarray, kpath: KPath, output: Path, options: PlotOptions,
               efermi: float = 0.0, label: str = "selection", notes: list[str] | None = None) -> Path:
    output = Path(output)
    options.validate(output)
    if options.emin is None or options.emax is None:
        raise ValueError("the energy window is not set; resolve it with resolve_window first")
    if energies.shape != weights.shape or energies.shape[0] != kpath.x.shape[0]:
        raise ValueError(f"energies {energies.shape}, weights {weights.shape} and {kpath.x.shape[0]} k-points do not match")
    shift = 0.0 if options.absolute else efermi
    energies = energies - shift
    fig = Figure(figsize=options.figsize)
    FigureCanvasAgg(fig)
    ax = fig.add_subplot(1, 1, 1)
    # Bands entirely below or entirely above the window cannot reach it; the rest form one contiguous slice.
    reach = np.flatnonzero((energies.max(axis=0) >= options.emin) & (energies.min(axis=0) <= options.emax))
    if reach.size:
        window = slice(int(reach[0]), int(reach[-1]) + 1)
        DRAW[options.mode](ax, fig, kpath, energies[:, window], weights[:, window], options, label)
    elif notes is not None:
        unit = "eV" if options.absolute else "eV relative to the Fermi level"
        notes.append(f"no bands in the window {options.emin:g} to {options.emax:g} {unit}; "
                     "the figure has empty axes and no legend")
    _decorate(ax, kpath, options, efermi - shift, legend=options.legend and bool(reach.size))
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=options.dpi)
    return output


# Axes limits, high-symmetry lines and labels, zero line at the Fermi level, title, and the legend when asked for.
def _decorate(ax: Axes, kpath: KPath, options: PlotOptions, zero: float, legend: bool = False) -> None:
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
    ax.set_ylabel("Energy (eV)" if options.absolute else r"$E - E_\mathrm{F}$ (eV)", fontsize=options.font_size)
    ax.tick_params(labelsize=options.font_size * 0.9, direction="in", top=False, right=True)
    ax.tick_params(axis="x", length=0)
    if options.title:
        ax.set_title(options.title, fontsize=options.font_size)
    if legend:
        ax.legend(loc="upper right", fontsize=options.font_size * 0.8, frameon=False, markerscale=0.6)
