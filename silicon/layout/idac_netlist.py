#!/usr/bin/env python3
"""FDVPD current DAC -- generate the reference netlist (LVS + simulation).

Architecture, from the paper's abstract: a 10 b differential current DAC with a
resistive load.  Each unit steers its current either into outp's load or into
a dump at VDD -- unipolar, not complementary into outn -- and a fixed block of
154 units (the "offset segments") always feeds outp:

    V_outn - V_outp = (c + N_OS) I_u R_D,        outn carries no current

outn exists only as vn's precharge path: an identical R_D from VDD, so both
sides of C_SAR see VDD through the same RC and supply noise cancels.  That
gives the detector exactly the model's lock point with nothing added: vn is
precharged to VDD and ramped down, so the residue is (c + N_OS) I_u R_D -
SR dt, zero at dt = T_frac + t_offset.  A complementary DAC would sit 400 mV
low (+-FS/2 about zero) and need an extend capacitor with a precision
reference, or a noisy offset sink, to put it back.

    segment   units    weight    count
    thermo      64      64 LSB     15      4 b MSB, thermometer
    binary   32..1   32..1 LSB      6      6 b LSB
                                  1023 units in the output array

    unit source   W/L = 0.84/24     390 nA, sigma(I)/I ~ 0.53 %  (idac_unit_size.py
                                    measured 0.75 % for 0.42/24 at half the current;
                                    same density, twice the area)
    switch pair   W/L = 0.42..1.68/0.15, 0.42 um per 16 units (switch_w)
                                    gate at VDD in saturation: it is also the
                                    source's cascode, so there is no separate one
    R_D           res_high_po_2p85, L = RD_L          ~2 kOhm, each side
    trim switch   W/L = 0.42/0.15   <= 12.5 uA: ~37 mV, negligible on a diode

The same-type reference diode is an array too: ``N = 64 + trim`` units, of which
``trim`` (7 b, 0..127) are switched in, so ``I_u = I_bias / N``.  That is the
R_D * C_SAR trim -- sky130 R*C spreads -30 %..+35 % (ngspice/rc_spread.sh) and
the background calibration only absorbs +-3 % -- done as a mirror ratio, which
scales every unit together and leaves linearity alone.

The reference is hierarchical, like the layout: idac instantiates idac_array,
which holds every unit source, diode and trim unit, and the dummies.

    python3 idac_netlist.py      # writes idac_ref.spice + idac_array_ref.spice
"""

from __future__ import annotations

import argparse
import pathlib

HERE = pathlib.Path(__file__).resolve().parent

NFET = "sky130_fd_pr__nfet_01v8"
RES = "sky130_fd_pr__res_high_po_2p85"

UNIT_W, UNIT_L = 0.84, 24
#: units that always feed outp: the model's V_OS margin, t_offset * SR =
#: 300 ps * 0.4 mV/ps = 120 mV = 154 x 0.78 mV.  Same unit devices, so the
#: offset tracks the DAC -- trim included.
N_OS = 154
OS_SEGMENTS = [64, 64, 26]
SW_W, SW_L = 0.42, 0.15
#: The load is high-sheet poly 2.85 um wide.  Its voltage coefficient is set
#: by the field along it, and the unipolar DAC puts up to 0.92 V across R_P
#: (the complementary one never exceeded 0.4 V, and its two sides' curvature
#: cancelled).  At 2 kOhm, INL came out 10.95 LSB with the 0.35 um resistor
#: (L = 1.24 um), 0.53 with 1.41 um (L = 7.6) and 0.077 with 2.85 um
#: (L = 16.4) -- against 0.035 with an ideal resistor.
RD_W = 2.85
RD_L = 16.43

THERMO_SEGMENTS = 15
THERMO_UNITS = 64
BINARY_UNITS = [32, 16, 8, 4, 2, 1]
REF_FIXED = 64
REF_TRIM = [1, 2, 4, 8, 16, 32, 64]     # 7 b

#: Per-segment cascode between the unit sources and the switch pair, placed
#: right at the switch.  The schematic did not need it (INL / DNL 0.069 / 0.040
#: LSB without, 0.062 / 0.008 with, at 20 ns per code) and it was first left
#: out.  Layout reversed that: each segment's source net is routed through
#: every channel column and across the array, 170 - 250 fF, and when the
#: segment switches sides that node moves ~20 mV -- charge the segment's own
#: current has to supply.  For the 1-unit segment that is ~17 ns, most of an
#: update period, and post-layout DNL went 0.04 -> 0.42 LSB.  With the cascode
#: next to the switch, the node that moves is small and local, and the big
#: routed one is the cascode's source, which holds still.
CASCODE = True
CASC_L = 0.5
CASC_NF = 2             # fingers per cascode
#: cascode width per unit.  Equal in every segment is what matters (same Vgs,
#: same unit Vds); 0.42 um at 390 nA runs it at twice the unit's density,
#: which cbias makes up for, and keeps a 64-unit finger at 13 um.
CASC_W_UNIT = 0.42

#: unit-array floorplan (layout/idac_array.py).  Here because the dummy count
#: it implies is part of the LVS reference.
ARRAY_COLS = 6
ARRAY_EDGE_DUMMIES = 2                  # stripes closing each column end


def segments() -> list[tuple[str, int]]:
    """The switched segments: 15 thermometer, 6 binary."""
    segs = [(f"t{k}", THERMO_UNITS) for k in range(THERMO_SEGMENTS)]
    segs += [(f"b{k}", n) for k, n in zip(range(5, -1, -1), BINARY_UNITS)]
    return segs


def switch_w(units: int) -> float:
    """Switch width for a segment: 0.42 um per 16 units, never below minimum.

    The switch node is 1.8 V - Vgs(switch) and it is the cascode's drain, so a
    switch run too hard starves the cascode: at 0.42 um for all, a thermometer
    switch carries 25 uA and at ss / -40 C left the cascode 49 mV of Vds.
    """
    return max(SW_W, round(SW_W * units / 16, 2))


def os_segments() -> list[tuple[str, int]]:
    """The offset units, in cascode-sized pieces: always on, into outp."""
    assert sum(OS_SEGMENTS) == N_OS
    return [(f"o{k}", n) for k, n in enumerate(OS_SEGMENTS)]


def array_units() -> list[tuple[str, int]]:
    """(net, units) for every unit in the array: outputs, diode, trim."""
    items = [(f"{'c' if CASCODE else 's'}_{n}", k)
             for n, k in segments() + os_segments()]
    items.append(("ibias", REF_FIXED))
    items += [(f"dt{k}", n) for k, n in enumerate(REF_TRIM)]
    return items


def array_rows() -> tuple[int, int]:
    """(inner rows, total rows) per column.  Inner rows hold every unit plus
    at least one spare pair, with an even total so every position has a
    point-symmetric partner."""
    n = sum(k for _, k in array_units())
    inner = -(-n // ARRAY_COLS)
    if inner * ARRAY_COLS == n:
        inner += 1
    if (inner * ARRAY_COLS) % 2:
        inner += 1
    return inner, inner + 2 * ARRAY_EDGE_DUMMIES


def array_dummies() -> int:
    return ARRAY_COLS * array_rows()[1] - sum(k for _, k in array_units())


def array_ports() -> list[str]:
    return ["vss"] + [n for n, _ in array_units()]


def cascode_w(units: int) -> float:
    """Per-finger cascode width: one unit width per unit of current, over
    CASC_NF fingers, never below the minimum.

    Equal current density in every cascode is what gives every unit source the
    same Vds.  A width clamped at the minimum for the small binary segments
    left the thermometer units at 0.31 V -- below V_dsat -- and the 1-unit
    segment's at 0.45 V: 0.41 LSB of static DNL, worse than no cascode.  Only
    the 1-unit segment is at the floor now (2x density), which moves one LSB
    by a few ppm.
    """
    return max(CASC_W_UNIT, round(CASC_W_UNIT * units / CASC_NF, 2))


def ports() -> list[str]:
    p = ["vdd", "vss", "ibias"] + (["cbias"] if CASCODE else []) + ["outp", "outn"]
    for name, _ in segments():
        p += [f"d_{name}", f"db_{name}"]
    p += [f"trim{k}" for k in range(len(REF_TRIM))]
    return p


def array_netlist() -> str:
    """The unit array alone: every source/diode unit, drains as ports."""
    L = [f".subckt idac_array {' '.join(array_ports())}"]
    for net, k in array_units():
        L.append(f"XA_{net} {net} ibias vss vss {NFET} W={UNIT_W} L={UNIT_L} m={k}")
    L.append("* dummies: drain and source to vss, gate on the shared sheet")
    L.append(f"XA_dum vss ibias vss vss {NFET} W={UNIT_W} L={UNIT_L} "
             f"m={array_dummies()}")
    L.append(".ends")
    return "\n".join(L) + "\n"


def netlist() -> str:
    L = []
    A = L.append
    A("* idac -- LVS reference for layout/idac.gds, generated by idac_netlist.py")
    A("* Do not edit by hand: change idac_netlist.py and regenerate.")
    A("")
    # hierarchical, like the layout: every unit source and diode unit lives in
    # idac_array, and idac adds the switches, trim switches and loads
    A(array_netlist())
    A(f".subckt idac {' '.join(ports())}")
    A("* the unit array: output segments, reference diode, trim units, dummies")
    A(f"Xarr {' '.join(array_ports())} idac_array")
    A("* load resistors: outp's carries the DAC current, outn's none -- it is")
    A("* vn's precharge path, matched so VDD noise cancels")
    A(f"XRp vdd outp vss {RES} L={RD_L}")
    A(f"XRn vdd outn vss {RES} L={RD_L}")
    A("* 7 b trim: each bit's units are diode-connected through a switch, so an")
    A("* unselected unit carries no current")
    for k in range(len(REF_TRIM)):
        A(f"XSt{k} ibias trim{k} dt{k} vss {NFET} W={SW_W} L={SW_L} m=1")
    A("* per segment: a switch pair steering it into outp or dumping it to VDD")
    for name, n in segments():
        src = f"s_{name}"
        if CASCODE:
            A(f"XC_{name} {src} cbias c_{name} vss {NFET} "
              f"W={cascode_w(n)} L={CASC_L} m={CASC_NF}")
        sw = switch_w(n)
        A(f"XP_{name} outp d_{name} {src} vss {NFET} W={sw} L={SW_L} m=1")
        A(f"XN_{name} vdd db_{name} {src} vss {NFET} W={sw} L={SW_L} m=1")
    A("* offset segments: always into outp, through a cascode like the rest")
    for name, n in os_segments():
        A(f"XC_{name} outp cbias c_{name} vss {NFET} "
          f"W={cascode_w(n)} L={CASC_L} m={CASC_NF}")
    A(".ends")
    return "\n".join(L) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", type=pathlib.Path, default=HERE / "idac_ref.spice")
    args = ap.parse_args()
    args.o.write_text(netlist(), encoding="utf-8", newline="\n")
    arr = args.o.with_name("idac_array_ref.spice")
    arr.write_text("* idac_array -- LVS reference, generated by idac_netlist.py\n"
                   + array_netlist(), encoding="utf-8", newline="\n")
    n_units = sum(n for _, n in segments())
    print(f"wrote {args.o.name}: {n_units} output units in {len(segments())} "
          f"segments, reference 64 + {sum(REF_TRIM)} trim, {len(ports())} ports")


if __name__ == "__main__":
    main()
