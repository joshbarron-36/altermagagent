"""Run the plotter as "python -m vasp_orbital_plotter" with the same arguments as vasp-orbital-bands."""

import sys

from .cli import main

sys.exit(main())
