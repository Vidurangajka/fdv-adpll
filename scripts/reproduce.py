"""Regenerate every figure, table and number in the report, then its PDF.

    python scripts/reproduce.py            # everything that runs on the host, ~15 min
    python scripts/reproduce.py --fast     # analytic figures, quick spec, no sims
    python scripts/reproduce.py --no-pdf   # skip printing docs/fdv-adpll-report.pdf

What each step regenerates, by where it appears in docs/report/report.html:

    tests       pytest (the fast subset with --fast)       Sec. 5 claims
    figures     scripts/run_all.py -> results/              Figs. 1-3, 5; Tables 1-2
    dataset     scripts/make_dataset.py --sim 48            Sec. 6, Fig. 4
    spec        silicon/spec/make_spec.py                   Tables 4-5
    report      scripts/build_report.py                     docs/fdv-adpll-report.pdf

The sky130A layout and circuit results (Sec. 7.3-7.4, Table 6, Fig. 7) need
the IIC-OSIC-TOOLS container and a LibreLane run, so they are not run from
here.  Their recorded outputs are committed (silicon/ngspice/*.json,
silicon/layout/*.gds) and silicon/README.md gives the command behind each
number; the checklist at the end of this script's output names them.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
PY = sys.executable

SILICON_STEPS = """\
Not run from here -- needs Docker and the IIC-OSIC-TOOLS image (silicon/osic.ps1
or silicon/osic.sh; see silicon/README.md, section by section):

  ramp sink     layout/verify_cell.sh ramp_sink ...; ngspice/corners.sh, ramp_se.sh
  current DAC   layout/idac_netlist.py + verify_cell.sh idac ...;
                ngspice/idac_tb.py, idac_corners.py, idac_mc.py -> idac_mc_loop.py
  bias          layout/bias.py; ngspice/bias_tb.py
  front end     layout/fe_netlist.py, sarfe.py; ngspice/fe_tb.py
  control       digital/ (LibreLane, config.json); tb_fdvpd_ctrl.v
  co-sim        ngspice/det_cosim.py
  top level     layout/verify_top.sh
"""


def step(name: str, cmd: list[str]) -> float:
    print(f"\n=== {name}: {' '.join(cmd)}", flush=True)
    t0 = time.time()
    env = dict(os.environ, MPLBACKEND="Agg")
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)
    dt = time.time() - t0
    print(f"=== {name} done in {dt:.0f} s", flush=True)
    return dt


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fast", action="store_true",
                    help="fast tests, analytic figures, --quick spec, no dataset sims")
    ap.add_argument("--no-pdf", action="store_true", help="do not rebuild the PDF")
    args = ap.parse_args()

    total = time.time()
    if args.fast:
        step("tests", [PY, "-m", "pytest", "-q", "-m", "not slow"])
        step("figures", [PY, "scripts/run_all.py", "--fast"])
        step("dataset", [PY, "scripts/make_dataset.py"])
        step("spec", [PY, "silicon/spec/make_spec.py", "--quick"])
    else:
        step("tests", [PY, "-m", "pytest", "-q"])
        step("figures", [PY, "scripts/run_all.py"])
        step("dataset", [PY, "scripts/make_dataset.py", "--sim", "48"])
        step("spec", [PY, "silicon/spec/make_spec.py"])
    if not args.no_pdf:
        step("report", [PY, "scripts/build_report.py"])

    print(f"\nall host steps done in {time.time() - total:.0f} s\n")
    print(SILICON_STEPS)


if __name__ == "__main__":
    main()
