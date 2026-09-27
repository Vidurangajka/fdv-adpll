#!/usr/bin/env python3
"""fdvpd_ctrl's pin template for LibreLane (FP_DEF_TEMPLATE).

The macro sits directly under the current DAC, its origin at the DAC's x = 0,
so every DAC switch and trim line is one straight m2 strip across the gap:

  north  trim[0..6], dsw[k] / dswb[k] -- each on the m2 track at or just left
         (dsw) / right (dswb) of the DAC's own pin, so a pair's two strips keep
         0.18 um between them where the DAC's pins are only 0.5 um apart
  east   the SAR front end's and the ramp's controls, upper half, m3
  south  everything that leaves the detector: REF, CKV, FCW, trim_in, the
         code, the frequency count

The DAC's pin x positions are read from its GDS (layout/build/idac/idac.gds),
so this has to run in the container:

    .\\osic.ps1 python3 digital/pins_def.py          # writes digital/pins.def
"""
from __future__ import annotations

import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
DAC_GDS = HERE.parent / "layout" / "build" / "idac" / "idac.gds"

DIE_W, DIE_H = 205.0, 110.0           # um
M2_P, M2_O = 0.46, 0.23               # sky130 hd met2 tracks (vertical)
M3_P, M3_O = 0.68, 0.34               # met3 tracks (horizontal)
PIN_W2, PIN_W3, PIN_L = 0.28, 0.60, 1.0   # m2 / m3 pin width, length into the die

#: dsw index -> DAC segment (dsw[20:15] = b0..b5, dsw[14:0] = t14..t0)
SEG = {k: f"t{k}" for k in range(15)} | {15 + j: f"b{5 - j}" for j in range(6)}

EAST = (["s1_n", "s2_n", "cmp_clk", "cmp_p", "cmp_n"]
        + [f"bpos[{k}]" for k in range(6)] + [f"bneg[{k}]" for k in range(6)]
        + [f"bmid[{k}]" for k in range(6)])
SOUTH = (["rst_n", "ref_clk", "ckv", "ckvd"] + [f"fcw_int[{k}]" for k in range(8)]
         + [f"fcw_frac[{k}]" for k in range(16)] + [f"trim_in[{k}]" for k in range(7)]
         + [f"code[{k}]" for k in range(7)] + ["code_valid", "railed"]
         + [f"count_ref[{k}]" for k in range(8)])
OUTPUTS = ({"ckvd", "s1_n", "s2_n", "cmp_clk", "code_valid", "railed"}
           | {f"{b}[{k}]" for b in ("bpos", "bneg", "bmid") for k in range(6)}
           | {f"code[{k}]" for k in range(7)} | {f"count_ref[{k}]" for k in range(8)}
           | {f"dsw[{k}]" for k in range(21)} | {f"dswb[{k}]" for k in range(21)}
           | {f"trim[{k}]" for k in range(7)})


def dac_pins() -> dict[str, float]:
    import klayout.db as kdb
    ly = kdb.Layout()
    ly.read(str(DAC_GDS))
    c = ly.cell("idac")
    out = {}
    for li in ly.layer_indexes():
        for s in c.shapes(li).each():
            if s.is_text():
                out[s.text_string] = s.text_pos.x / 1000.0
    return out


def track(x: float, pitch: float, off: float, side: str) -> float:
    n = (x - off) / pitch
    n = math.floor(n + 1e-9) if side == "left" else math.ceil(n - 1e-9)
    return round(off + n * pitch, 3)


def north() -> list[tuple[str, float]]:
    dp = dac_pins()
    pins = [(f"trim[{k}]", track(dp[f"trim{k}"], M2_P, M2_O, "left")) for k in range(7)]
    for k in range(21):
        s = SEG[k]
        pins.append((f"dsw[{k}]", track(dp[f"d_{s}"], M2_P, M2_O, "left")))
        pins.append((f"dswb[{k}]", track(dp[f"db_{s}"], M2_P, M2_O, "right")))
    xs = sorted(x for _, x in pins)
    assert all(b - a >= M2_P - 1e-6 for a, b in zip(xs, xs[1:])), "north pins collide"
    return pins


def um(v: float) -> int:
    return int(round(v * 1000))


def rows() -> list[tuple[str, str, int, int, str]]:
    """Every pin as (name, layer, x, y, side), nm, the macro's own coordinates.

    layout/top.py routes to these, so it and the macro cannot disagree.
    """
    out = []
    # every pin lies inside the die, flush with its edge, as OpenROAD's own do
    for name, x in north():
        out.append((name, "met2", um(x), um(DIE_H - PIN_L / 2), "N"))
    # east: from 1 um under the top, every other m3 track
    y = track(DIE_H - 3.0, M3_P, M3_O, "left")
    for name in EAST:
        out.append((name, "met3", um(DIE_W - PIN_L / 2), um(y), "E"))
        y = round(y - 2 * M3_P, 3)
    # south: spread over the width, on m2 tracks
    step = (DIE_W - 10.0) / (len(SOUTH) - 1)
    for i, name in enumerate(SOUTH):
        out.append((name, "met2", um(track(5.0 + i * step, M2_P, M2_O, "left")),
                    um(PIN_L / 2), "S"))
    return out


def main() -> None:
    rows_ = rows()
    L = ["VERSION 5.8 ;", 'DIVIDERCHAR "/" ;', 'BUSBITCHARS "[]" ;',
         "DESIGN fdvpd_ctrl ;", "UNITS DISTANCE MICRONS 1000 ;",
         f"DIEAREA ( 0 0 ) ( {um(DIE_W)} {um(DIE_H)} ) ;", f"PINS {len(rows_)} ;"]
    hl = um(PIN_L / 2)
    for name, layer, x, y, side in rows_:
        d = "OUTPUT" if name in OUTPUTS else "INPUT"
        # m3 pins at m3's 0.3 um minimum width would be refused: 0.6, as OpenROAD
        hw = um((PIN_W2 if layer == "met2" else PIN_W3) / 2)
        box = f"( -{hw} -{hl} ) ( {hw} {hl} )" if side in "NS" else f"( -{hl} -{hw} ) ( {hl} {hw} )"
        L += [f"    - {name} + NET {name} + DIRECTION {d} + USE SIGNAL",
              f"      + PORT", f"        + LAYER {layer} {box}",
              f"        + FIXED ( {x} {y} ) N ;"]
    L += ["END PINS", "END DESIGN", ""]
    out = HERE / "pins.def"
    out.write_text("\n".join(L), newline="\n")
    print(f"wrote {out.name}: {len(rows_)} pins, die {DIE_W} x {DIE_H} um")


if __name__ == "__main__":
    main()
