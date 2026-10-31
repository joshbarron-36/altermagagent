"""Stage 3 of the orbital plotter: spin panels, signed magnetisation blocks, overlays, labels, styles and large-file options."""

import dataclasses
import json
import os
import shutil
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "fixtures"))

import orbital_vasprun_builder as fx  # noqa: E402
from matplotlib.collections import PathCollection  # noqa: E402
from matplotlib.colors import to_hex, to_rgba  # noqa: E402
from vasp_orbital_plotter import plotting, projections  # noqa: E402
from vasp_orbital_plotter.kpath import (KPath, build_kpath, parse_label_map, rename_labels, stride_indices,  # noqa: E402
                                        subsample, without_labels)
from vasp_orbital_plotter import cli  # noqa: E402
from vasp_orbital_plotter.kpath import split_label_map  # noqa: E402
from vasp_orbital_plotter.options import PANEL_SIZE, KPathError, OptionsError, PlotOptions, load_rcparams  # noqa: E402
from vasp_orbital_plotter.parser import parse_vasprun  # noqa: E402
from vasp_orbital_plotter.plotting import (NEGATIVE_COLOUR, NEUTRAL_GREY, POSITIVE_COLOUR, SERIES_COLOURS, Panel,  # noqa: E402
                                           Series, grey_centred_cmap)
from vasp_orbital_plotter.projections import (Block, SelectionError, bands_reaching, build_selection,  # noqa: E402
                                              build_series, choose_blocks, select_weights, state_totals)

from test_orbital_plotting import _Base, png_size  # noqa: E402
from test_orbital_projections import hand_sum  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "orbital"
REAL_RUN = Path("/home/clawc/development/osi-hackathon/data/diana_MnTe/hackaton_exa/9.2-band/vasprun.xml")
TOL = 1e-9


# Parse one synthetic fixture by its short name.
def fixture(name: str):
    return parse_vasprun(FIXTURES / f"{name}_vasprun.xml")


# A three-point k-path on one branch with the given tick labels at x = 0, 1, 2.
def small_path(labels: list[str]) -> KPath:
    return KPath(x=np.array([0.0, 1.0, 2.0]), branches=[(0, 2)], segment_starts=[],
                 ticks=[(float(i), lab) for i, lab in enumerate(labels)], segment_source="test", label_source="test")


class BlockChoice(unittest.TestCase):
    """--spin and --component pick the right projection and eigenvalue blocks, and refuse layouts they do not fit."""

    # Parse the three layouts once.
    @classmethod
    def setUpClass(cls) -> None:
        cls.nonmag, cls.collinear, cls.soc = fixture("nonmagnetic"), fixture("collinear_f"), fixture("soc")

    # Nonmagnetic files have one block; --spin and a magnetisation component are refused.
    def test_nonmagnetic(self) -> None:
        choice = choose_blocks(self.nonmag)
        self.assertEqual((choice.blocks, choice.signed, choice.caption), ([Block(0, None)], False, None))
        for kwargs, text in (({"spin": "up"}, "nonmagnetic"), ({"component": "mz"}, "only block is total")):
            with self.subTest(kwargs=kwargs), self.assertRaises(SelectionError) as caught:
                choose_blocks(self.nonmag, **kwargs)
            self.assertIn(text, str(caught.exception))

    # Collinear files draw both spins by default, each with its own eigenvalue block; one spin on request.
    def test_collinear(self) -> None:
        both = choose_blocks(self.collinear)
        self.assertEqual(both.blocks, [Block(0, "spin up"), Block(1, "spin down")])
        self.assertIn("two panels", both.note)
        down = choose_blocks(self.collinear, spin="down")
        self.assertEqual((down.blocks, down.eigen(down.blocks[0]), down.tag), ([Block(1, None)], 1, "down"))
        self.assertEqual(both.tag, "both")
        self.assertIsNone(choose_blocks(self.collinear, spin="up").note)
        with self.assertRaises(SelectionError) as caught:
            choose_blocks(self.collinear, component="mx")
        self.assertIn("use --spin", str(caught.exception))

    # Noncollinear files take --component; mx, my, mz are signed with the unverified-sign caption; --spin is refused.
    def test_noncollinear(self) -> None:
        self.assertEqual(choose_blocks(self.soc).blocks, [Block(0, None)])
        self.assertEqual((choose_blocks(self.soc).tag, choose_blocks(self.soc).quantity), ("", "projection weight"))
        self.assertIsNone(choose_blocks(self.soc).caption)
        for index, name in enumerate(("mx", "my", "mz"), start=1):
            choice = choose_blocks(self.soc, component=name)
            self.assertEqual((choice.blocks, choice.signed, choice.label), ([Block(index, None)], True, name))
            self.assertEqual((choice.eigen(choice.blocks[0]), choice.tag, choice.quantity), (0, name, name))
            self.assertIn("sign convention unverified", choice.caption)
            self.assertIn("LSORBIT = T", choice.caption)
        with self.assertRaises(SelectionError) as caught:
            choose_blocks(self.soc, spin="up")
        self.assertIn("use --component", str(caught.exception))
        for kwargs, text in (({"normalize": True}, "drop --normalize"), ({"overlay": True}, "--overlay is not available")):
            with self.subTest(kwargs=kwargs), self.assertRaises(SelectionError) as caught:
                choose_blocks(self.soc, component="mz", **kwargs)
            self.assertIn(text, str(caught.exception))
        with self.assertRaises(SelectionError) as caught:
            choose_blocks(self.soc, component="mw")
        self.assertIn("must be one of", str(caught.exception))

    # A noncollinear run without SOC is plotted like an SOC run, and its caption states LSORBIT = F.
    def test_lsorbit_false(self) -> None:
        nosoc = dataclasses.replace(self.soc, parameters={**self.soc.parameters, "LSORBIT": False})
        self.assertIn("LSORBIT = F", choose_blocks(nosoc).caption)
        self.assertIn("LSORBIT = F", choose_blocks(nosoc, component="mz").caption)


class SpinWeights(_Base):
    """Collinear weights and eigenvalues come from the chosen spin block; both spins make two panels."""

    # Parse the collinear fixture.
    def setUp(self) -> None:
        super().setUp()
        self.data = fixture("collinear_f")

    # Run the command line with plot_panels replaced by a recorder; return the recorded arguments.
    def captured(self, *args: str) -> tuple[tuple, dict]:
        with patch.object(plotting, "plot_panels", return_value=self.out / "x.png") as mock:
            code, _, err = self.run_cli("--vasprun", str(FIXTURES / "collinear_f_vasprun.xml"), *args,
                                        "--emin", "-9", "--emax", "9", "--output", str(self.out / "x.png"))
        self.assertEqual(code, 0, err)
        return mock.call_args.args, mock.call_args.kwargs

    # --spin down uses projection block 1 and eigenvalue block 1, which differ from spin up in the fixture.
    def test_spin_down(self) -> None:
        self.assertGreater(np.abs(self.data.eigenvalues[0] - self.data.eigenvalues[1]).max(), 0.1)
        (panels, *_), _ = self.captured("--atoms", "3", "--orbitals", "f", "--spin", "down")
        self.assertEqual(len(panels), 1)
        np.testing.assert_allclose(panels[0].energies, self.data.eigenvalues[1] - self.data.efermi, atol=TOL)
        np.testing.assert_allclose(panels[0].series[0].weights, hand_sum([2], list(range(9, 16)), c=1), atol=TOL)

    # --spin both gives two panels, up then down, each with its own energies and weights.
    def test_spin_both(self) -> None:
        (panels, *_), _ = self.captured("--atoms", "1", "--orbitals", "s", "--spin", "both")
        self.assertEqual([p.heading for p in panels], ["spin up", "spin down"])
        for c, panel in enumerate(panels):
            np.testing.assert_allclose(panel.energies, self.data.eigenvalues[c] - self.data.efermi, atol=TOL)
            np.testing.assert_allclose(panel.series[0].weights, hand_sum([0], [0], c=c), atol=TOL)

    # Without --figsize the figure is one default panel width per panel; an explicit --figsize is kept as given.
    def test_two_panel_size(self) -> None:
        for extra, size in (((), (240, 90)), (("--figsize", "6,4.5"), (120, 90))):
            with self.subTest(extra=extra):
                out = self.out / f"both{len(extra)}.png"
                code, text, err = self.run_cli("--vasprun", str(FIXTURES / "collinear_f_vasprun.xml"), "--orbitals", "d",
                                               "--mode", "color", "--dpi", "20", *extra, "--output", str(out))
                self.assertEqual(code, 0, err)
                self.assertEqual(png_size(out), size)
                self.assertIn("note: collinear ISPIN = 2 file: spin up and spin down are drawn as two panels", text)

    # Both panels share the energy axis; only the first shows the energy label.
    def test_shared_axis(self) -> None:
        kpath = build_kpath(self.data)
        sel = build_selection(self.data, orbitals="d")
        panels = [Panel(self.data.eigenvalues[c], [Series(select_weights(self.data, sel, c), "d")], f"spin {s}")
                  for c, s in enumerate(("up", "down"))]
        fig = plotting._figure(panels, kpath, PlotOptions(emin=-3, emax=3, figsize=(12.0, 4.5)), 0.0, [], False,
                               "projection weight", None)
        first, second = fig.axes[:2]
        self.assertIs(second.get_shared_y_axes().joined(first, second), True)
        self.assertEqual(second.get_ylabel(), "")
        self.assertEqual([first.get_title(), second.get_title()], ["spin up", "spin down"])


class SignedComponents(_Base):
    """mx, my and mz keep their sign through the sum and are drawn with signed colours and a caption."""

    # Synthetic SOC fixture: mz of Cr2 d is the signed hand sum of block 3.
    def test_soc_fixture_weights(self) -> None:
        data = fixture("soc")
        sel = build_selection(data, atoms="2", orbitals="d")
        for c in (1, 2, 3):
            np.testing.assert_allclose(select_weights(data, sel, c), hand_sum([1], [4, 5, 6, 7, 8], c=c, signed=True),
                                       atol=TOL)
        self.assertTrue((select_weights(data, sel, 3) < 0).all())

    # Trimmed real MnTe fixture: the block sums match direct sums of the parsed array, signs included.
    def test_trimmed_mnte(self) -> None:
        data = fixture("mnte_9.2-band_trimmed")
        sel = build_selection(data, elements="Mn", orbitals="d")
        for c in range(4):
            expected = data.projections[c][:, :, 0:2, 4:9].sum(axis=(2, 3))
            np.testing.assert_allclose(select_weights(data, sel, c), expected, atol=TOL)
        state = data.projections.sum(axis=(3, 4))
        self.assertTrue((np.sqrt((state[1:] ** 2).sum(axis=0)) <= state[0] + 1e-3).all())

    # Colour mode: symmetric norm about zero, grey middle, colour bar named after the component; the caption is printed.
    def test_color_cli(self) -> None:
        out = self.out / "mz.png"
        code, text, err = self.run_cli("--vasprun", str(FIXTURES / "soc_vasprun.xml"), "--orbitals", "d", "--component", "mz",
                                       "--mode", "color", "--dpi", "30", "--emin", "-9", "--emax", "9", "--output", str(out))
        self.assertEqual(code, 0, err)
        self.assertIn("block mz", text)
        self.assertIn("note: mz block as read from vasprun.xml, sign convention unverified; LSORBIT = T", text)
        self.assertTrue(out.is_file())
        np.testing.assert_allclose(grey_centred_cmap("RdBu_r")(0.5), to_rgba(NEUTRAL_GREY), atol=0.01)

    # The norm of signed colour mode is centred on zero and symmetric; limits on one side of zero are refused.
    def test_signed_norm(self) -> None:
        panel = Panel(np.array([[0.0, 1.0]]), [Series(np.array([[-0.2, 0.5]]), "x")])
        norm = plotting._color_norm([panel], PlotOptions(emin=-3, emax=3), signed=True)
        self.assertEqual((norm.vmin, norm.vcenter, norm.vmax), (-0.5, 0.0, 0.5))
        with self.assertRaises(OptionsError):
            plotting._color_norm([panel], PlotOptions(emin=-3, emax=3, cmin=0.1, cmax=0.5), signed=True)
        self.assert_user_error("--vasprun", str(FIXTURES / "soc_vasprun.xml"), "--component", "mz", "--mode", "color",
                               "--cmin", "0.1", "--cmax", "0.5", "--output", str(self.out / "x.png"),
                               contains="both sides of zero")

    # Fat-band mode: red markers for positive values, blue for negative, area proportional to |value|, legend per sign.
    def test_signed_fatband(self) -> None:
        energies = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
        weights = np.array([[0.5, -0.25], [0.5, -0.25], [0.0, 0.0]])
        options = PlotOptions(emin=-3, emax=3, marker_scale=100.0, figsize=PANEL_SIZE)
        fig = plotting._figure([Panel(energies, [Series(weights, "Mn d")])], small_path(["Γ", "C₀", "Y"]), options, 0.0, [],
                               True, "mz", "caption text")
        ax = fig.axes[0]
        scatters = [c for c in ax.collections if isinstance(c, PathCollection) and c.get_offsets().size]
        colours = {to_hex(c.get_facecolor()[0]): sorted(set(c.get_sizes())) for c in scatters}
        self.assertEqual(colours, {POSITIVE_COLOUR: [50.0], NEGATIVE_COLOUR: [25.0]})
        self.assertEqual([t.get_text() for t in ax.get_legend().get_texts()], ["Mn d: mz > 0", "Mn d: mz < 0"])
        self.assertEqual(fig.get_supxlabel(), "caption text")


@unittest.skipUnless(REAL_RUN.is_file(), "Diana's MnTe 9.2-band run is not on this machine")
class RealBlocks(_Base):
    """Sanity of the four noncollinear blocks of the real MnTe 9.2-band run, and a signed plot from it."""

    # Parse the real file once.
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = parse_vasprun(REAL_RUN, progress=False)

    # The total per state, summed over ions and orbitals, lies in [0, 1.05], and |m| never exceeds it.
    def test_block_sanity(self) -> None:
        state = self.data.projections.sum(axis=(3, 4))
        total, m = state[0], np.sqrt((state[1:] ** 2).sum(axis=0))
        self.assertGreaterEqual(total.min(), 0.0)
        self.assertLessEqual(total.max(), 1.05)
        self.assertTrue((m <= total + 1e-3).all(), float((m - total).max()))

    # Summed over occupied states, Mn1's mx, my point along its MAGMOM and Mn2's the opposite way (consistency, not proof).
    def test_moment_along_magmom(self) -> None:
        magmom = np.array(self.data.parameters["MAGMOM"]).reshape(-1, 3)
        occupied = self.data.occupations[0][:, :, None]
        for ion in (0, 1):
            moment = np.array([(self.data.projections[c][:, :, ion, :].sum(axis=-1) * occupied[..., 0]).sum()
                               for c in (1, 2, 3)])
            cosine = moment @ magmom[ion] / (np.linalg.norm(moment) * np.linalg.norm(magmom[ion]))
            self.assertGreater(cosine, 0.99, (ion, moment))

    # Signed mz colour plot of the real run, written to a temporary folder only.
    def test_real_mz_plot(self) -> None:
        out = self.out / "mz.png"
        code, text, err = self.run_cli("--vasprun", str(REAL_RUN), "--element", "Mn", "--orbitals", "d", "--component", "mz",
                                       "--mode", "color", "--dpi", "50", "--output", str(out))
        self.assertEqual(code, 0, err)
        self.assertIn("sign convention unverified", text)
        self.assertEqual(png_size(out), (300, 225))


class Overlays(_Base):
    """--overlay draws each --select term as its own series with its own colour and legend entry."""

    # build_series: one selection per term, names from --label-series, and its refusals.
    def test_build_series(self) -> None:
        data = fixture("soc")
        series = build_series(data, ["Cr:d", "Se:p"])
        self.assertEqual([s.label for s in series], ["Cr d", "Se p"])
        self.assertEqual([s.label for s in build_series(data, ["Cr:d", "Se:p"], ["metal", "anion"])], ["metal", "anion"])
        for args, text in (((None,), "needs --select"), ((["Cr:d"], ["a", "b"]), "one name per term"),
                           ((["Cr:s"] * 7,), "at most 6")):
            with self.subTest(text=text), self.assertRaises(SelectionError) as caught:
                build_series(data, *args)
            self.assertIn(text, str(caught.exception))

    # Fat-band overlay: one colour per series in the fixed order and one legend entry each.
    def test_fatband_overlay(self) -> None:
        with patch.object(plotting, "plot_panels", return_value=self.out / "x.png") as mock:
            code, text, err = self.run_cli("--vasprun", str(FIXTURES / "soc_vasprun.xml"), "--select", "Cr:d", "--select",
                                           "Se:p", "--overlay", "--label-series", "Cr d", "--label-series", "Se p",
                                           "--output", str(self.out / "x.png"))
        self.assertEqual(code, 0, err)
        panel = mock.call_args.args[0][0]
        self.assertEqual([(s.label, s.colour) for s in panel.series], [("Cr d", SERIES_COLOURS[0]), ("Se p", SERIES_COLOURS[1])])
        self.assertIn("selection: Cr d;", text)
        self.assertIn("selection: Se p;", text)
        fig = plotting._figure([panel], build_kpath(fixture("soc")),
                               PlotOptions(emin=-9, emax=9, figsize=PANEL_SIZE), 1.23456789, [], False, "projection weight", None)
        self.assertEqual([t.get_text() for t in fig.axes[0].get_legend().get_texts()], ["Cr d", "Se p"])

    # Colour-mode overlay: each state takes its dominant series' colour, full at the top weight, grey at zero.
    def test_dominant_colour(self) -> None:
        energies = np.array([[0.0, 1.0], [0.0, 1.0]])
        first = Series(np.array([[1.0, 0.0], [1.0, 0.0]]), "A", SERIES_COLOURS[0])
        second = Series(np.array([[0.2, 0.0], [0.2, 0.0]]), "B", SERIES_COLOURS[1])
        two_points = KPath(x=np.array([0.0, 1.0]), branches=[(0, 1)], segment_starts=[], ticks=[(0.0, ""), (1.0, "")],
                           segment_source="t", label_source="t")
        fig = plotting._figure([Panel(energies, [first, second])], two_points,
                               PlotOptions(mode="color", emin=-3, emax=3, figsize=PANEL_SIZE), 0.0, [], False, "projection weight",
                               None)
        lines = next(c for c in fig.axes[0].collections if len(c.get_segments()) == 4)
        colours = {to_hex(c) for c in lines.get_colors()}
        self.assertEqual(colours, {SERIES_COLOURS[0], NEUTRAL_GREY})
        self.assertEqual(len(fig.axes), 1)
        self.assertEqual([t.get_text() for t in fig.axes[0].get_legend().get_texts()], ["A", "B"])

    # Overlay plots write files in both modes on the real-layout fixture.
    def test_overlay_files(self) -> None:
        for mode in ("fatband", "color"):
            with self.subTest(mode=mode):
                out = self.out / f"overlay_{mode}.png"
                code, _, err = self.run_cli("--vasprun", str(FIXTURES / "mnte_9.2-band_trimmed_vasprun.xml"), "--select",
                                            "Mn:d", "--select", "Te:p", "--overlay", "--mode", mode, "--emin", "-60",
                                            "--emax", "0", "--dpi", "30", "--output", str(out))
                self.assertEqual(code, 0, err)
                self.assertTrue(out.is_file())


class Labels(_Base):
    """--label-map renames path labels with mathtext, --no-labels removes them, and mistakes are refused."""

    # Keys match a label as printed or as typed; parts of joined labels are renamed one by one.
    def test_rename(self) -> None:
        kpath = rename_labels(small_path(["Γ", "C₀", "Y|Σ₀"]), parse_label_map("G=A,C0=C$_0$,Σ₀=S"))
        self.assertEqual([lab for _, lab in kpath.ticks], ["A", "C$_0$", "Y|S"])
        self.assertIn("renamed by --label-map", kpath.label_source)
        self.assertEqual([lab for _, lab in without_labels(small_path(["Γ", "M", "K"])).ticks], ["", "", ""])

    # Unknown keys, malformed entries, odd $ signs and bad mathtext raise KPathError.
    def test_label_errors(self) -> None:
        with self.assertRaises(KPathError) as caught:
            rename_labels(small_path(["Γ", "M", "K"]), parse_label_map("X=Z"))
        self.assertIn("labels on the path: Γ, M, K", str(caught.exception))
        for text in ("G", "=A", "G=C$_0", r"G=$\notacommand{x}$"):
            with self.subTest(text=text), self.assertRaises(KPathError):
                parse_label_map(text)

    # Renamed and subscript labels render in a saved figure; --no-labels prints its label source.
    def test_cli(self) -> None:
        fixture_path = str(FIXTURES / "soc_vasprun.xml")
        code, text, err = self.run_cli("--vasprun", fixture_path, "--label-map", "G=Γ$_1$,M=M$_0$", "--dpi", "30",
                                       "--output", str(self.out / "renamed.svg"))
        self.assertEqual(code, 0, err)
        self.assertIn("labels from vasprun <kpoints_labels>, renamed by --label-map", text)
        code, text, err = self.run_cli("--vasprun", fixture_path, "--no-labels", "--dpi", "30",
                                       "--output", str(self.out / "none.png"))
        self.assertEqual(code, 0, err)
        self.assertIn("labels from none (--no-labels)", text)


class Formatting(_Base):
    """Styles, the rcParams file, frame lines, inward ticks and vector output."""

    # The talk style scales font size, line width and marker area; paper changes nothing.
    def test_styles(self) -> None:
        self.assertEqual(PlotOptions().styled(), PlotOptions())
        talk = PlotOptions(style="talk", font_size=10.0).styled()
        self.assertEqual((talk.font_size, talk.linewidth, talk.marker_scale), (14.0, 0.8 * 1.5, 120.0))

    # The rcParams file reaches the figure; invalid JSON, a non-object and unknown keys are refused.
    def test_rcparams_file(self) -> None:
        good = self.out / "rc.json"
        good.write_text(json.dumps({"axes.linewidth": 3.25}))
        self.assertEqual(load_rcparams(good), {"axes.linewidth": 3.25})
        out = self.out / "rc.svg"
        code, _, err = self.run_cli("--vasprun", str(FIXTURES / "soc_vasprun.xml"), "--rcparams-file", str(good),
                                    "--output", str(out))
        self.assertEqual(code, 0, err)
        self.assertIn("stroke-width: 3.25", out.read_text())
        for content, text in (("{", "not valid JSON"), ("[1]", "one JSON object"), ('{"axes.linewidthx": 1}', "axes.linewidthx")):
            bad = self.out / "bad.json"
            bad.write_text(content)
            with self.subTest(content=content):
                self.assert_user_error("--vasprun", str(FIXTURES / "soc_vasprun.xml"), "--rcparams-file", str(bad),
                                       "--output", str(self.out / "x.png"), contains=text)

    # --despine hides the top and right frame lines; ticks point inwards on the band axes and on the colour bar.
    def test_despine_and_ticks(self) -> None:
        data = fixture("soc")
        weights = select_weights(data, build_selection(data, orbitals="d"), 0)
        options = PlotOptions(mode="color", emin=-9, emax=9, despine=True, figsize=PANEL_SIZE)
        fig = plotting._figure([Panel(data.eigenvalues[0], [Series(weights, "d")])], build_kpath(data), options, 1.2, [],
                               False, "projection weight", None)
        ax, bar = fig.axes
        self.assertEqual((ax.spines["top"].get_visible(), ax.spines["right"].get_visible()), (False, False))
        self.assertEqual(ax.yaxis.get_major_ticks()[0]._tickdir, "in")
        self.assertEqual(bar.yaxis.get_major_ticks()[0]._tickdir, "in")

    # SVG and PDF are vector by default; --rasterize embeds an image.
    def test_vector_output(self) -> None:
        for suffix, marker in ((".svg", "<image"), (".pdf", "/Subtype /Image")):
            for raster in (False, True):
                with self.subTest(suffix=suffix, raster=raster):
                    out = self.out / f"v{int(raster)}{suffix}"
                    code, _, err = self.run_cli("--vasprun", str(FIXTURES / "soc_vasprun.xml"), "--orbitals", "d",
                                                "--emin", "-9", "--emax", "9", "--dpi", "40", *(("--rasterize",) if raster else ()),
                                                "--output", str(out))
                    self.assertEqual(code, 0, err)
                    self.assertEqual(marker in out.read_bytes().decode("latin-1"), raster)


class LargeFiles(_Base):
    """--kpoint-stride keeps segment ends, the band window cuts the weight sum, and --float32 halves the memory."""

    # Every stride-th point of each segment plus each segment's end; branches and starts are renumbered.
    def test_stride_indices(self) -> None:
        self.assertEqual(stride_indices(10, [5], 2).tolist(), [0, 2, 4, 5, 7, 9])
        self.assertEqual(stride_indices(10, [5], 1).tolist(), list(range(10)))
        with self.assertRaises(KPathError):
            stride_indices(10, [], 0)
        kpath = KPath(x=np.arange(10.0), branches=[(0, 4), (5, 9)], segment_starts=[5], ticks=[(0.0, "A")],
                      segment_source="t", label_source="t")
        thin = subsample(kpath, stride_indices(10, [5], 2))
        self.assertEqual((thin.x.tolist(), thin.branches, thin.segment_starts), ([0, 2, 4, 5, 7, 9], [(0, 2), (3, 5)], [3]))

    # Weights over a band window and kept k-points equal the full computation on those states.
    def test_window_and_kpoints(self) -> None:
        data = fixture("collinear_f")
        sel = build_selection(data, elements="Cr", orbitals="d")
        full = select_weights(data, sel, 1, normalize=True)
        part = select_weights(data, sel, 1, normalize=True, kpoints=np.array([0, 3]), bands=slice(1, 3))
        np.testing.assert_allclose(part, full[[0, 3], 1:3], atol=TOL)
        self.assertEqual(bands_reaching([np.array([[-5.0, 0.0, 5.0], [-4.0, 1.0, 6.0]])], -1.0, 1.0), slice(1, 2))
        self.assertEqual(bands_reaching([np.array([[-5.0, 9.0]])], -1.0, 1.0), slice(0, 0))

    # The stride and the band window cut the array the einsum sees, and the normalising total can be passed in.
    def test_stride_reaches_einsum(self) -> None:
        data = fixture("collinear_f")
        sel = build_selection(data, atoms="1", orbitals="d")
        seen: list[tuple[int, ...]] = []
        einsum = np.einsum

        # Record the shape of the projection operand, then sum as numpy would.
        def recording(spec: str, *operands: np.ndarray) -> np.ndarray:
            seen.append(operands[0].shape)
            return einsum(spec, *operands)

        kept = np.array([0, 3])
        with patch.object(projections.np, "einsum", recording):
            full = select_weights(data, sel, 0)
            thin = select_weights(data, sel, 0, kpoints=kept, bands=slice(1, 3))
        self.assertEqual([shape[:3] for shape in seen], [(data.nkpoints, data.nbands, 1), (2, 2, 1)])
        self.assertEqual(thin.shape, (2, 2))
        uneven = np.array([0, 1, 3])
        np.testing.assert_allclose(select_weights(data, sel, 0, kpoints=uneven), full[uneven], atol=TOL)
        np.testing.assert_allclose(thin, full[[0, 3], 1:3], atol=TOL)
        total = state_totals(data, 0, kept, slice(1, 3))
        np.testing.assert_allclose(select_weights(data, sel, 0, normalize=True, kpoints=kept, bands=slice(1, 3), total=total),
                                   thin / total, atol=TOL)

    # The command line reports the kept k-points and the bands that reach the window.
    def test_cli_reports(self) -> None:
        code, text, err = self.run_cli("--vasprun", str(FIXTURES / "soc_vasprun.xml"), "--kpoint-stride", "2",
                                       "--emin", "-9", "--emax", "9", "--dpi", "30", "--output", str(self.out / "s.png"))
        self.assertEqual(code, 0, err)
        self.assertIn("k-points: 4 of 4 kept (--kpoint-stride 2, segment ends kept)", text)
        self.assertIn(f"bands: {fx.NBANDS} of {fx.NBANDS} reach the window", text)

    # --float32 stores single-precision projections equal to the double-precision ones within the file's 4 decimals.
    def test_float32(self) -> None:
        single = parse_vasprun(FIXTURES / "soc_vasprun.xml", float32=True)
        double = fixture("soc")
        self.assertEqual(single.projections.dtype, np.float32)
        self.assertEqual(single.projections.nbytes * 2, double.projections.nbytes)
        np.testing.assert_allclose(single.projections, double.projections, atol=5e-5)
        code, _, err = self.run_cli("--vasprun", str(FIXTURES / "soc_vasprun.xml"), "--float32", "--dpi", "30",
                                    "--output", str(self.out / "f32.png"))
        self.assertEqual(code, 0, err)


class Refusals(_Base):
    """Every new refusal is one "error:" line with exit code 2."""

    # Each refused combination names the problem.
    def test_refusals(self) -> None:
        nonmag, collinear = str(FIXTURES / "nonmagnetic_vasprun.xml"), str(FIXTURES / "collinear_f_vasprun.xml")
        soc = str(FIXTURES / "soc_vasprun.xml")
        out = ("--output", str(self.out / "x.png"))
        cases = [
            ((nonmag, "--spin", "up"), "nonmagnetic"),
            ((soc, "--spin", "down"), "use --component"),
            ((nonmag, "--component", "mx"), "only block is total"),
            ((collinear, "--component", "mz"), "use --spin"),
            ((soc, "--spin", "sideways"), "invalid choice"),
            ((soc, "--component", "mw"), "invalid choice"),
            ((soc, "--component", "mz", "--normalize"), "drop --normalize"),
            ((soc, "--component", "mz", "--overlay", "--select", "Cr:d"), "--overlay is not available"),
            ((soc, "--overlay"), "needs --select"),
            ((soc, "--overlay", "--element", "Cr", "--select", "Cr:d"), "with --overlay"),
            ((soc, "--label-series", "A", "--select", "Cr:d"), "add --overlay"),
            ((soc, "--overlay", "--select", "Cr:d", "--label-series", "A", "--label-series", "B"), "one name per term"),
            ((soc, "--overlay", *sum((("--select", "Cr:s") for _ in range(7)), ())), "at most 6"),
            ((soc, "--label-map", "Q=Z"), "not on the path"),
            ((soc, "--label-map", "G=$_"), "odd number of $"),
            ((soc, "--label-map", "G=A", "--no-labels"), "use one of them"),
            ((soc, "--kpoint-stride", "0"), "1 or more"),
            ((soc, "--style", "poster"), "invalid choice"),
            ((soc, "--cmap-signed", "NoSuchMap"), "not known to Matplotlib"),
            ((soc, "--rcparams-file", str(self.out / "missing.json")), "cannot be read"),
        ]
        for args, text in cases:
            with self.subTest(args=args):
                self.assert_user_error("--vasprun", *args, *out, contains=text)


class ReviewFixes(_Base):
    """Fixes of the stage 3 code review: label-map splitting, pre-parse refusals, signed legends and panel checks."""

    # Commas inside $...$ mathtext or inside braces do not split a --label-map entry.
    def test_label_map_commas(self) -> None:
        self.assertEqual(parse_label_map("X=X$_{1,2}$,G=Γ"), [("X", "X$_{1,2}$"), ("G", "Γ")])
        self.assertEqual(split_label_map("A=a{1,2}, B=b ,"), ["A=a{1,2}", "B=b"])

    # Run the command line with parse_vasprun replaced by a mock; return the error text and whether it was called.
    def refused_before_parse(self, *args: str, contains: str) -> bool:
        with patch.object(cli, "parse_vasprun") as parse:
            self.assert_user_error(*args, contains=contains)
        return parse.called

    # In the vasprun's own folder without --output, the run is refused before the file is parsed and nothing is written.
    def test_default_output_before_parse(self) -> None:
        shutil.copy(FIXTURES / "nonmagnetic_vasprun.xml", self.out / "vasprun.xml")
        cwd = os.getcwd()
        self.addCleanup(os.chdir, cwd)
        os.chdir(self.out)
        self.assertFalse(self.refused_before_parse("--orbitals", "d", contains="pass --output"))
        self.assertEqual(sorted(p.name for p in self.out.iterdir()), ["vasprun.xml"])

    # Refusals that need no file come before the parse: stride, overlay without terms, series limit, label count.
    def test_flag_refusals_before_parse(self) -> None:
        soc, out = str(FIXTURES / "soc_vasprun.xml"), ("--output", str(self.out / "x.png"))
        cases = [
            (("--kpoint-stride", "0"), "--kpoint-stride must be 1 or more, not 0"),
            (("--overlay",), "--overlay needs --select"),
            (("--overlay", *sum((("--select", "Cr:s") for _ in range(7)), ())), "at most 6 series, but 7"),
            (("--overlay", "--select", "Cr:d", "--label-series", "A", "--label-series", "B"), "2 --label-series names for 1"),
        ]
        for args, text in cases:
            with self.subTest(args=args):
                self.assertFalse(self.refused_before_parse("--vasprun", soc, *args, *out, contains=text))
        self.assertFalse((self.out / "x.png").exists())

    # Signed colour mode draws no legend entry, even with --legend; the colour bar names the quantity.
    def test_signed_color_legend(self) -> None:
        panel = Panel(np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]), [Series(np.array([[0.5, -0.2]] * 3), "Mn d")])
        options = PlotOptions(mode="color", emin=-3, emax=3, legend=True, figsize=PANEL_SIZE)
        self.assertEqual(plotting._legend_handles(panel, options, True, "mz"), [])
        fig = plotting._figure([panel], small_path(["Γ", "M", "K"]), options, 0.0, [], True, "mz", None)
        self.assertIsNone(fig.axes[0].get_legend())
        self.assertEqual(fig.axes[1].get_ylabel(), "mz projection (sign convention unverified)")

    # A panel without series and a panel of the wrong k-point count are refused before drawing; nothing is written.
    def test_panel_checks(self) -> None:
        path, options, out = small_path(["Γ", "M", "K"]), PlotOptions(emin=-3, emax=3), self.out / "p.png"
        good = Panel(np.zeros((3, 2)), [Series(np.ones((3, 2)), "a")])
        cases = [([good, Panel(np.zeros((3, 2)), [])], "panel 2 has no series"),
                 ([good, Panel(np.zeros((2, 2)), [Series(np.ones((2, 2)), "b")])], "panel 2 has energies for 2 k-points"),
                 ([Panel(np.zeros((3, 2)), [Series(np.ones((3, 1)), "c")])], "of 'c' do not match"), ([], "no panels")]
        for panels, text in cases:
            with self.subTest(text=text), self.assertRaises(ValueError) as caught:
                plotting.plot_panels(panels, path, out, options)
            self.assertIn(text, str(caught.exception))
        self.assertFalse(out.exists())
        with self.assertRaises(ValueError):
            plotting._figure([good], path, options, 0.0, [], False, "projection weight", None)


if __name__ == "__main__":
    unittest.main()
