#!/usr/bin/env python3
"""Take the DAC Monte-Carlo through the loop: fractional-spur yield.

    python silicon/ngspice/idac_mc_loop.py          # after idac_mc.py, ~30 min

For every draw in idac_mc.json (21 segment currents from sky130 Monte-Carlo),
the behavioural model's DAC is rebuilt with exactly those segment weights and
the sky130 operating point is run in the spec's fractional channel (bit 5).
The weights are normalised so full scale is 1023 LSB: gain is the trim's job,
and what is left is the DAC's non-linearity -- the thing the tolerance budget's
three-seed dac_sigma_lsb number was a stand-in for.

Also reports each draw's static INL / DNL, and one run with ideal weights for
the floor the rest of the loop sets.

Writes ngspice/idac_mc_loop.json.
"""

from __future__ import annotations

import json
import pathlib
import statistics as st
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE.parent / "spec"))

SPUR_TARGET = -60.0
N_CYCLES = 1 << 16


def weights(draw: dict, units: dict) -> tuple[list, list]:
    """(15 thermometer weights, 6 binary weights LSB first), in LSB, full
    scale normalised to 1023."""
    seg = draw["segments"]
    unit = sum(seg.values()) / 1023.0
    thermo = [seg[f"t{k}"] / unit for k in range(15)]
    binary = [seg[f"b{k}"] / unit for k in range(6)]
    return thermo, binary


def levels(thermo, binary) -> np.ndarray:
    codes = np.arange(1024)
    cum = np.concatenate(([0.0], np.cumsum(thermo)))
    bits = ((codes & 63)[:, None] >> np.arange(6)) & 1
    return cum[codes >> 6] + bits @ np.asarray(binary)


def inl_dnl(lv) -> tuple[float, float]:
    gain = (lv[-1] - lv[0]) / 1023
    fit = lv[0] + gain * np.arange(1024)
    return (float(np.max(np.abs(lv - fit)) / gain),
            float(np.max(np.abs(np.diff(lv) / gain - 1))))


def run(args):
    idx, thermo, binary = args
    import fdvadpll.blocks as blocks
    from fdvadpll import FdvPll
    from make_spec import sky130_design, fractional, DISCARD_FRACTION

    # the original class, pinned once per worker: each call re-patches, and
    # subclassing the previous call's patch would inherit its weights
    base = blocks.__dict__.setdefault("_OrigCurrentDAC", blocks.CurrentDAC)

    class MeasuredDAC(base):
        def __init__(self, params, rng):
            super().__init__(params, rng)          # sigma 0: no draw
            if thermo is not None:
                self.w_thermo = np.asarray(thermo)
                self.w_bin = np.asarray(binary)
                self.levels = levels(thermo, binary) * (1.0 + params.dac_gain_err)
                self.frac = self.levels / self.n_codes

    blocks.CurrentDAC = MeasuredDAC
    d = fractional(sky130_design(), 5)
    res = FdvPll(d, seed=idx).run(N_CYCLES, discard=N_CYCLES // DISCARD_FRACTION)
    sat = float(res.saturated[res.discard:].mean())
    locked = bool(sat < 1e-3 and res.cycle_slips == 0)
    return idx, float(res.fractional_spur()[1]), float(res.worst_spur()[1]), locked


def main():
    mc = json.loads((HERE / "idac_mc.json").read_text())
    draws = mc["draws"]
    jobs = [(i, *weights(d, mc["units"])) for i, d in enumerate(draws)]
    lin = [inl_dnl(levels(t, b)) for _, t, b in jobs]
    with ProcessPoolExecutor(4) as ex:
        ideal = list(ex.map(run, [(10_000 + s, None, None) for s in range(4)]))
        out = list(ex.map(run, jobs))
    spurs = [o[1] for o in out]
    locked = [o[3] for o in out]
    good = sum(1 for s, l in zip(spurs, locked) if l and s <= SPUR_TARGET)
    ideal_spur = max(o[1] for o in ideal)
    print(f"draws                 {len(out)}")
    print(f"ideal-DAC floor       {ideal_spur:.1f} dBc (worst of {len(ideal)} seeds)")
    print(f"fractional spur, bit 5: median {st.median(spurs):.1f}, "
          f"95th {np.percentile(spurs, 95):.1f}, worst {max(spurs):.1f} dBc")
    print(f"locked                {sum(locked)}/{len(out)}")
    print(f"yield vs {SPUR_TARGET:.0f} dBc     {good}/{len(out)} = {100*good/len(out):.1f} %")
    inl = [x[0] for x in lin]
    dnl = [x[1] for x in lin]
    print(f"static INL (LSB)      median {st.median(inl):.3f}, 95th {np.percentile(inl, 95):.3f}")
    print(f"static DNL (LSB)      median {st.median(dnl):.3f}, 95th {np.percentile(dnl, 95):.3f}")
    (HERE / "idac_mc_loop.json").write_text(json.dumps({
        "target_dbc": SPUR_TARGET, "n_cycles": N_CYCLES, "ideal_floor_dbc": ideal_spur,
        "draws": [{"spur_dbc": o[1], "worst_spur_dbc": o[2], "locked": o[3],
                   "inl": l[0], "dnl": l[1]} for o, l in zip(out, lin)]}, indent=1))


if __name__ == "__main__":
    main()
