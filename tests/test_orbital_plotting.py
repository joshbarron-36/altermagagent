"""Stage 2 of the orbital plotter: figure files from both modes, and the command line with its exit codes."""

import contextlib
import io
import os
import shutil
import struct
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

from vasp_orbital_plotter import cli
from vasp_orbital_plotter.options import OptionsError, PlotOptions

from _fixtures import REPO, tmp_folder

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "orbital"
REAL_RUN = Path("/home/clawc/development/osi-hackathon/data/diana_MnTe/hackaton_exa/9.2-band/vasprun.xml")
MIN_BYTES = 5_000


# Width and height in pixels from the IHDR chunk of a PNG file.
def png_size(path: Path) -> tuple[int, int]:
    head = path.read_bytes()[:24]
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        raise AssertionError(f"{path} is not a PNG file")
    return struct.unpack(">II", head[16:24])


class _Base(unittest.TestCase):
    """Runs the command line with captured output inside a fresh temporary folder."""

    # Make a temporary output folder that is removed after the test.
    def setUp(self) -> None:
        self.out = tmp_folder(self, "orbital_plot_")

    # Run the command line; return exit code, stdout and stderr.
    def run_cli(self, *args: str) -> tuple[int, str, str]:
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            code = cli.main(list(args))
        return code, stdout.getvalue(), stderr.getvalue()

    # Assert a user error: exit code 2, one "error:" line on stderr, no traceback.
    def assert_user_error(self, *args: str, contains: str = "") -> str:
        code, _, err = self.run_cli(*args)
        self.assertEqual(code, 2, err)
        lines = err.strip().splitlines()
        self.assertEqual(len(lines), 1, err)
        self.assertTrue(lines[0].startswith("error: "), err)
        self.assertNotIn("Traceback", err)
        self.assertIn(contains, err)
        return err


class PlotFiles(_Base):
    """Both modes write PNG, SVG and PDF files of the expected kind and size."""

    # Fat-band PNG of the nonmagnetic fixture: 4 x 3 inches at 50 dpi gives 200 x 150 pixels.
    def test_fatband_png(self) -> None:
        out = self.out / "fat.png"
        code, text, err = self.run_cli("--vasprun", str(FIXTURES / "nonmagnetic_vasprun.xml"), "--atoms", "1,3",
                                       "--orbitals", "d", "--figsize", "4,3", "--dpi", "50", "--output", str(out))
        self.assertEqual(code, 0, err)
        self.assertIn(f"wrote {out}", text)
        self.assertGreater(out.stat().st_size, 2_000)
        self.assertEqual(png_size(out), (200, 150))

    # Colour mode with a colour bar, at the default 300 dpi and 6 x 4.5 inches.
    def test_color_png(self) -> None:
        out = self.out / "colour.png"
        code, _, err = self.run_cli("plot", "--vasprun", str(FIXTURES / "soc_vasprun.xml"), "--element", "Cr",
                                    "--orbitals", "t2g", "--mode", "color", "--output", str(out))
        self.assertEqual(code, 0, err)
        self.assertGreater(out.stat().st_size, MIN_BYTES)
        self.assertEqual(png_size(out), (1800, 1350))

    # SVG and PDF are chosen by the extension.
    def test_svg_and_pdf(self) -> None:
        for name, magic in (("bands.svg", b"<svg"), ("bands.pdf", b"%PDF")):
            with self.subTest(name=name):
                out = self.out / name
                code, _, err = self.run_cli("--vasprun", str(FIXTURES / "collinear_f_vasprun.xml"), "--select", "Cr:f",
                                            "--select", "2:p", "--mode", "color", "--normalize", "--output", str(out))
                self.assertEqual(code, 0, err)
                self.assertGreater(out.stat().st_size, MIN_BYTES)
                self.assertIn(magic, out.read_bytes()[:400])

    # Absolute energies, a manual Fermi energy, labels and a legend are accepted together.
    def test_options(self) -> None:
        out = self.out / "abs.png"
        code, text, err = self.run_cli("--vasprun", str(FIXTURES / "nonmagnetic_vasprun.xml"), "--absolute", "--efermi", "0.5",
                                       "--high-symmetry", "1:G,4:K", "--legend", "--no-title", "--marker-scale", "200",
                                       "--alpha", "0.3", "--dpi", "40", "--output", str(out))
        self.assertEqual(code, 0, err)
        self.assertIn("labels from --high-symmetry option", text)
        self.assertGreater(png_size(out)[0], 0)

    # Collinear files print the spin-up note.
    def test_collinear_note(self) -> None:
        code, text, _ = self.run_cli("--vasprun", str(FIXTURES / "collinear_f_vasprun.xml"), "--orbitals", "s",
                                     "--dpi", "30", "--output", str(self.out / "c.png"))
        self.assertEqual(code, 0)
        self.assertIn("note: collinear ISPIN = 2 file", text)

    # A window that no band reaches still writes the figure, exits 0, prints a note and draws no legend.
    def test_empty_window(self) -> None:
        out = self.out / "empty.png"
        code, text, err = self.run_cli("--vasprun", str(FIXTURES / "nonmagnetic_vasprun.xml"), "--emin", "50", "--emax", "60",
                                       "--legend", "--dpi", "30", "--output", str(out))
        self.assertEqual(code, 0, err)
        self.assertEqual(err, "")
        self.assertIn("note: no bands in the window 50 to 60 eV relative to the Fermi level", text)
        self.assertTrue(out.is_file())

    # Diana's spelling "--atom 2 --orbital pz" gives the same selection as "--atoms 2 --orbitals pz".
    def test_singular_spellings(self) -> None:
        texts = []
        for atoms, orbitals in (("--atom", "--orbital"), ("--atoms", "--orbitals")):
            with self.subTest(atoms=atoms):
                code, text, err = self.run_cli("--vasprun", str(FIXTURES / "nonmagnetic_vasprun.xml"), atoms, "2", orbitals,
                                               "pz", "--dpi", "30", "--output", str(self.out / f"{atoms[2:]}.png"))
                self.assertEqual(code, 0, err)
                texts.append(next(line for line in text.splitlines() if line.startswith("selection:")))
        self.assertEqual(texts[0], texts[1])
        self.assertIn("selection: Se pz; 1 ions, 1 orbitals", texts[0])


class CommandLine(_Base):
    """Exit codes, listings, defaults and the rule never to write next to the vasprun.xml by default."""

    # --list-atoms and --list-orbitals work as flags of the main command and of inspect.
    def test_listings(self) -> None:
        fixture = str(FIXTURES / "nonmagnetic_vasprun.xml")
        code, text, _ = self.run_cli("--vasprun", fixture, "--list-atoms")
        self.assertEqual(code, 0)
        self.assertIn("     2  Se", text)
        code, text, _ = self.run_cli("--vasprun", fixture, "--list-orbitals")
        self.assertEqual(code, 0)
        self.assertIn("     9  x2-y2", text)
        self.assertIn("eg = dz2+x2-y2", text)
        code, text, _ = self.run_cli("inspect", "--vasprun", fixture, "--list-orbitals", "--list-atoms")
        self.assertEqual(code, 0)
        self.assertIn("number  element", text)
        self.assertIn("t2g = dxy+dxz+dyz", text)
        code, text, _ = self.run_cli("inspect", "--vasprun", fixture)
        self.assertEqual(code, 0)
        self.assertIn("layout: nonmagnetic", text)

    # User errors give exit code 2 and one "error:" line.
    def test_user_errors(self) -> None:
        fixture = str(FIXTURES / "nonmagnetic_vasprun.xml")
        out = str(self.out / "x.png")
        self.assert_user_error("--vasprun", fixture, "--orbitals", "dq", "--output", out, contains="available orbitals: s, py")
        self.assert_user_error("--vasprun", fixture, "--atoms", "0", "--output", out, contains="VASP numbering starts at 1")
        self.assert_user_error("--vasprun", fixture, "--element", "Fe", "--output", out, contains="elements present: Cr, Se")
        self.assert_user_error("--vasprun", str(self.out / "missing.xml"), contains="no file at")
        self.assert_user_error("--vasprun", fixture, "--emin", "2", "--emax", "1", "--output", out, contains="lower energy limit 2")
        self.assert_user_error("--vasprun", fixture, "--output", str(self.out / "x.jpg"), contains=".png")
        self.assert_user_error("--vasprun", fixture, "--mode", "pie", contains="mode")
        self.assert_user_error("--vasprun", fixture, "--no-such-option", contains="unrecognized")
        self.assert_user_error("--vasprun", fixture, "--high-symmetry", "0:G", "--output", out, contains="numbering starts at 1")
        self.assert_user_error("--vasprun", fixture, "--cmap", "nonsense", "--output", out, contains="colour map")
        self.assert_user_error("--vasprun", str(FIXTURES / "no_projected_vasprun.xml"), contains="LORBIT = 11")
        self.assert_user_error("--vasprun", str(FIXTURES / "truncated_vasprun.xml"), contains="well-formed")
        self.assert_user_error("--vasprun", fixture, "--mode", "colour", contains="invalid choice")
        self.assert_user_error("--vasprun", fixture, "--font-size", "0", "--output", out, contains="font size must be positive")
        self.assert_user_error("--vasprun", fixture, "--cmin", "1", "--cmax", "1", "--output", out, contains="colour limit 1")
        empty = self.out / "empty_field.xml"
        head, tail = (FIXTURES / "nonmagnetic_vasprun.xml").read_text(encoding="latin-1").split("<projected>")
        empty.write_text(head + "<projected>" + tail.replace("<field>   py</field>", "<field></field>", 1), encoding="latin-1")
        self.assert_user_error("--vasprun", str(empty), "--orbitals", "p", "--output", out, contains="empty <field>")
        self.assert_user_error("--vasprun", str(empty), "--list-orbitals", contains="empty <field>")
        for alias in ("--ato", "--elements", "--orbit", "--normalise", "--labels", "-o"):
            with self.subTest(alias=alias):
                self.assert_user_error("--vasprun", fixture, alias, "1", "--output", out, contains="unrecognized")
        self.assert_user_error("inspect", fixture, contains="unrecognized")

    # --list-atoms needs no <projected> block; --list-orbitals on such a file is a one-line user error.
    def test_listings_without_projections(self) -> None:
        fixture = str(FIXTURES / "no_projected_vasprun.xml")
        for prefix in ((), ("inspect",)):
            with self.subTest(prefix=prefix):
                code, text, err = self.run_cli(*prefix, "--vasprun", fixture, "--list-atoms")
                self.assertEqual(code, 0, err)
                self.assertIn("number  element", text)
                self.assert_user_error(*prefix, "--vasprun", fixture, "--list-orbitals", contains="no orbital labels")

    # Only the project's own errors become "error:" lines; any other ValueError or OSError keeps its traceback.
    def test_unexpected_errors_propagate(self) -> None:
        fixture = str(FIXTURES / "nonmagnetic_vasprun.xml")
        for error in (ValueError("internal"), OSError("disk")):
            with self.subTest(error=type(error).__name__), patch.object(cli, "build_selection", side_effect=error):
                with self.assertRaises(type(error)):
                    self.run_cli("--vasprun", fixture, "--output", str(self.out / "x.png"))
        self.assert_user_error("--vasprun", str(self.out), "--list-atoms", contains="")

    # PlotOptions.validate refuses a font size of 0 or less, a negative line width and colour limits out of order.
    def test_validate(self) -> None:
        PlotOptions(linewidth=0.0, cmin=0.0, cmax=0.5).validate()
        for bad, text in ((PlotOptions(font_size=0), "font size"), (PlotOptions(font_size=-2), "font size"),
                          (PlotOptions(linewidth=-0.1), "line width"), (PlotOptions(cmin=0.5, cmax=0.5), "colour limit"),
                          (PlotOptions(cmin=0.6, cmax=0.5), "colour limit")):
            with self.subTest(text=text), self.assertRaises(OptionsError) as caught:
                bad.validate()
            self.assertIn(text, str(caught.exception))

    # Every PlotOptions field but the parsed rcparams has a flag of the same name whose default is the PlotOptions default.
    def test_defaults_from_plot_options(self) -> None:
        import dataclasses
        from vasp_orbital_plotter.options import PlotOptions
        self.assertEqual(cli.options_from_args(cli.plot_parser().parse_args([])), PlotOptions())
        flagged = [f for f in dataclasses.fields(PlotOptions) if f.metadata.get("flag", True)]
        self.assertEqual([f.name for f in dataclasses.fields(PlotOptions) if f not in flagged], ["rcparams"])
        self.assertTrue(all(hasattr(cli.plot_parser().parse_args([]), f.name) for f in flagged))

    # Listings and inspect never import Matplotlib, the plotting module or the k-path module.
    def test_listings_import_no_matplotlib(self) -> None:
        code = ("import sys; from vasp_orbital_plotter import cli; "
                f"cli.main(['--vasprun', {str(FIXTURES / 'nonmagnetic_vasprun.xml')!r}, '--list-atoms', '--list-orbitals']); "
                f"cli.main(['inspect', '--vasprun', {str(FIXTURES / 'soc_vasprun.xml')!r}]); "
                "print(sorted(m for m in sys.modules if m.startswith(('matplotlib', 'vasp_orbital_plotter.plotting', "
                "'vasp_orbital_plotter.kpath'))))")
        result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True,
                                env={**os.environ, "PYTHONPATH": str(REPO)})
        self.assertTrue(result.stdout.strip().endswith("[]"), result.stdout + result.stderr)

    # Run in a folder with vasprun.xml: --vasprun defaults to it, but the default output name is refused there.
    def test_default_vasprun_and_refusal(self) -> None:
        shutil.copy(FIXTURES / "nonmagnetic_vasprun.xml", self.out / "vasprun.xml")
        cwd = os.getcwd()
        self.addCleanup(os.chdir, cwd)
        os.chdir(self.out)
        self.assert_user_error("--atoms", "1", "--orbitals", "pz", contains="pass --output")
        self.assertEqual(sorted(p.name for p in self.out.iterdir()), ["vasprun.xml"])
        code, text, err = self.run_cli("--atoms", "1", "--orbitals", "pz", "--dpi", "30", "--output", "atom1_pz.pdf")
        self.assertEqual(code, 0, err)
        self.assertTrue((self.out / "atom1_pz.pdf").is_file())

    # Elsewhere, the default output is a name built from the selection in the current folder.
    def test_default_name(self) -> None:
        cwd = os.getcwd()
        self.addCleanup(os.chdir, cwd)
        os.chdir(self.out)
        code, text, err = self.run_cli("--vasprun", str(FIXTURES / "nonmagnetic_vasprun.xml"), "--element", "Cr",
                                       "--orbitals", "d", "--dpi", "30")
        self.assertEqual(code, 0, err)
        self.assertTrue((self.out / "bands_Cr_d_fatband.png").is_file(), text)


@unittest.skipUnless(REAL_RUN.is_file(), "Diana's MnTe 9.2-band run is not on this machine")
class RealFigures(_Base):
    """The real MnTe run plots in both modes; files are written to a temporary folder only."""

    # Mn d fat bands and Te p colour, with segment and label sources reported.
    def test_real_run(self) -> None:
        for mode, element, orbitals in (("fatband", "Mn", "d"), ("color", "Te", "p")):
            with self.subTest(mode=mode):
                out = self.out / f"{element}_{orbitals}_{mode}.png"
                code, text, err = self.run_cli("--vasprun", str(REAL_RUN), "--element", element, "--orbitals", orbitals,
                                               "--mode", mode, "--dpi", "100", "--output", str(out))
                self.assertEqual(code, 0, err)
                self.assertIn("k-path segments from KPOINTS; labels from KLABELS", text)
                self.assertEqual(png_size(out), (600, 450))
                self.assertGreater(out.stat().st_size, 20_000)


if __name__ == "__main__":
    unittest.main()
