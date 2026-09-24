"""Run every figure script in order, into ``results/``.

    python scripts/run_all.py            # everything
    python scripts/run_all.py --fast     # skip the long simulation sweeps

The analytic-only scripts take a second or two; the ones that run the
event-driven simulator take a few minutes.
"""

from __future__ import annotations

import argparse
import importlib
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

#: (module, needs a full simulation sweep)
SCRIPTS = [
    ("design_summary", False),
    ("fig05_noise_vs_power", False),
    ("fig12_noise_breakdown", False),
    ("fig17_power_and_fom", False),
    ("compare_architectures", False),
    ("fig11_phase_noise", True),
    ("calibration_demo", True),
    ("fig15_spur_vs_fraction", True),
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fast", action="store_true",
                    help="skip the scripts that run the event-driven simulator")
    args = ap.parse_args()

    total = time.time()
    for name, slow in SCRIPTS:
        if slow and args.fast:
            print(f"\n[skip] {name}  (simulation sweep)")
            continue
        print(f"\n[{name}]")
        t0 = time.time()
        importlib.import_module(name).main()
        print(f"  ({time.time() - t0:.1f} s)")
    print(f"\nall done in {time.time() - total:.1f} s -> results/")


if __name__ == "__main__":
    main()
