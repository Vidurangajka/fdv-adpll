#!/usr/bin/env python3
"""FDVPD front end -- generate the reference netlist of cell ``sarfe``.

Everything between the DAC, the ramp sink and the control block:

    S1p, S1n   encode switches: DAC outp -> vp, outn -> vn   (PMOS on s1_n;
               vp's is a transmission gate, its NMOS on a local inverter)
    Sr, Sd     the ramp: the sink's current steered into vn (Sr) or VDD (Sd),
               NMOS, gates from vg-powered drivers on s2_n -- see steer()
    CDAC       per side: a fixed MIM and six switched weights, 0.5 pF in all
    comparator StrongARM, NMOS input, with inverting output buffers

The residue is vn - vp.  The CDAC's weights 32..1 LSB (LSB = 195 uV) are MIM
units of the PDK's minimum size -- 2 x 2 um, ~9.3 fF -- switched about a mid
level by ladder taps rather than built from tiny capacitors:

    weight   units   step about vrm      vd step = 2 * n * C_u * dV / C_tot
      32       8      +-21.0 mV
      16       4      +-21.0
       8       2      +-21.0
       4       1      +-21.0
       2       1      +-10.5
       1       1      +- 5.25

Each bottom plate is three-state: mid (vrm), H or L.  bpos[k] puts the p plate
at H and the n plate at L, which lowers vn - vp; bneg[k] the opposite; bmid[k]
(neither, decoded by the control block) puts both at mid.  The first comparison takes the sign with
every plate at mid, so the 7 b result needs only six switched weights.

The comparator's outputs are inverted: both low in reset, and cmp_p rises when
vn > vp.

    python3 fe_netlist.py        # writes sarfe_ref.spice
"""
from __future__ import annotations

import argparse
import os
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
NFET = "sky130_fd_pr__nfet_01v8"
PFET = "sky130_fd_pr__pfet_01v8"
MIM = "sky130_fd_pr__cap_mim_m3_1"

C_TOT = 0.5e-12
UNIT_W = UNIT_L = 2.0              # um: the MIM minimum, ~9.3 fF
#: sky130 MIM, typical: 2.00 fF/um^2 and 0.19 fF/um of perimeter, on a plate
#: the model shrinks by m3_dw = 25 nm a side


def c_mim(w: float, l: float) -> float:
    w, l = w - 0.025, l - 0.025
    return 2.00e-15 * w * l + 0.19e-15 * 2 * (w + l)


C_UNIT = c_mim(UNIT_W, UNIT_L)
#: (units, H tap, L tap) per weight, LSB first
WEIGHTS = [(1, "vrp5", "vrn5"), (1, "vrp10", "vrn10"), (1, "vrp20", "vrn20"),
           (2, "vrp20", "vrn20"), (4, "vrp20", "vrn20"), (8, "vrp20", "vrn20")]
#: parasitics expected on each top plate: switches, comparator input, wiring
C_PAR = 15e-15
#: the fixed MIM, square, sized so the side totals C_TOT
_c_fix = C_TOT - C_PAR - sum(n for n, _, _ in WEIGHTS) * C_UNIT
# 2 s^2 + 0.76 s = C in fF, s the drawn side less m3_dw
FIX_SIDE = round((-0.76 + (0.76 ** 2 + 8 * _c_fix * 1e15) ** 0.5) / 4 + 0.025, 2)


#: steering switch: 10 um, two fingers
SW_W, SW_M = 5, 2
#: tracking bias for its gate: a 1/10 replica over an 18 um / 0.5 um diode
REP_W, DIODE_W, DIODE_M, DIODE_L = 1, 4.5, 4, 0.5
#: vg's decap, an NMOS gate: 2 x 20 um x 5 um, ~1.8 pF
CAP_W, CAP_L, CAP_M = 20, 5, 2


def steer() -> list[str]:
    """The ramp switch and its gate bias (device lines, subcircuit nets).

    A PMOS in series with the ramp starves the sink once vn falls past ~1.1 V
    (vn ends at 0.88 V on the longest window); an NMOS with its gate at VDD
    drops into triode near the end, and its channel charge then lands on vn
    (ramp_nl2 ~0.005 /V against a ~0.0035 limit).  So Sr's gate sits lower,
    at vg: the replica stack puts Sr's source ~0.65 V above ground whatever
    the corner, which leaves the sink its headroom and keeps Sr saturated
    down to vn ~0.8 V.  nl2 0.0020 - 0.0031 /V over ss/tt/ff x -40..85 C.

    The drivers make before break at the window's start (Sr on, then Sd
    off); at its end Sd comes back one gate later, which only moves rs.
    vg is also the bias current's input: 20 uA from the bias block.
    """
    return [
        f"XSr vn g rs vss {NFET} W={SW_W} L=0.15 m={SW_M}",
        f"XSd vdd gb rs vss {NFET} W={SW_W} L=0.15 m={SW_M}",
        f"XDp1 g s2_n vg vdd {PFET} W=2 L=0.15",
        f"XDn1 g s2_n vss vss {NFET} W=0.5 L=0.15",
        f"XDp2 gb g vg vdd {PFET} W=2 L=0.15",
        f"XDn2 gb g vss vss {NFET} W=0.5 L=0.15",
        f"XRs vg vg m vss {NFET} W={REP_W} L=0.15",
        f"XRd m m vss vss {NFET} W={DIODE_W} L={DIODE_L} m={DIODE_M}",
        f"XCg vss vg vss vss {NFET} W={CAP_W} L={CAP_L} m={CAP_M}",
    ]


#: encode switches, um: PMOS on both sides, and an NMOS beside vp's only
#: (FE_S1_WN overrides the NMOS for sizing sweeps)
S1_W = 5
S1_WN = float(os.environ.get("FE_S1_WN", 5))
#: the NMOS half is low-Vt: near vp = 1 V both halves have ~0.2 V of
#: overdrive, and more width buys conductance only with more channel charge
#: dumped on vp at turn-off (FE_S1_NMOS overrides it for sweeps)
NFET_LVT = "sky130_fd_pr__nfet_01v8_lvt"
S1_NMOS = os.environ.get("FE_S1_NMOS", NFET_LVT)


def encode() -> list[str]:
    """The encode switches: outp -> vp and outn -> vn, closed while s1_n is low.

    vp has to follow the DAC down to 0.88 V, where a lone PMOS has |Vgs|
    0.9 V against a threshold its body effect (the nwell at VDD) has raised to
    about the same: the full-range front-end test (fe_tb.py --start 0 --step
    121) read a 3.9 mV rms bow off it, against 0.12 mV with ideal switches.
    So vp's switch is a transmission gate, its NMOS half carrying the bottom
    of the range.  Same test, residual rms: 0.59 mV with a 5 um standard NMOS
    (worst at vp ~ 1.0 V, where neither half has overdrive), 0.49 at 10 um
    (the extra channel charge now errs at the bottom), 0.23 mV with a 5 um
    low-Vt NMOS and 0.33 at 10 um -- 5 um low-Vt it is.

    vn's is not.  vn is only ever precharged to VDD, which the PMOS passes
    well, and an NMOS there would do nothing but hang a nonlinear junction on
    the node the ramp sweeps: with transmission gates on both sides the bow
    grew with their width (0.44 mV rms at 10 um, 1.09 at 20) and the ramp
    slowed.  vp holds still while the ramp runs.  s1 comes from an inverter
    here, not a second pin.
    """
    return [
        f"XS1p outp s1_n vp vdd {PFET} W={S1_W} L=0.15",
        f"XS1n outn s1_n vn vdd {PFET} W={S1_W} L=0.15",
        f"XS1pn outp s1 vp vss {S1_NMOS} W={S1_WN:g} L=0.15",
        f"XIs1p s1 s1_n vdd vdd {PFET} W=2 L=0.15",
        f"XIs1n s1 s1_n vss vss {NFET} W=1 L=0.15",
    ]


def ports() -> list[str]:
    p = ["vdd", "vss", "outp", "outn", "rs", "vg", "s1_n", "s2_n", "cmp_clk",
         "cmp_p", "cmp_n", "vp", "vn"]
    p += [f"bpos{k}" for k in range(6)] + [f"bneg{k}" for k in range(6)]
    p += [f"bmid{k}" for k in range(6)]
    p += ["vrp20", "vrp10", "vrp5", "vrm", "vrn5", "vrn10", "vrn20"]
    return p


def netlist() -> str:
    L = []
    A = L.append
    A("* sarfe -- LVS reference for layout/sarfe.gds, generated by fe_netlist.py")
    A(f".subckt sarfe {' '.join(ports())}")
    A("* encode switches: transmission gates")
    L.extend(encode())
    A("* the ramp switch: current steering, gate bias tracking the corner")
    L.extend(steer())
    A("* CDAC: fixed part, top plate to the side, bottom to vss")
    for side in ("p", "n"):
        A(f"XCf{side} v{side} vss {MIM} W={FIX_SIDE} L={FIX_SIDE}")
    # one instance per MIM unit: the model's mf scales only its mismatch, not
    # its capacitance, and m= would give the units one shared mismatch draw
    A("* switched weights: MIM units, three-state bottom plates")
    for k, (n, th, tl) in enumerate(WEIGHTS):
        sw_w = 0.84 if n <= 2 else 1.68
        A(f"* weight {2**k}: {n} unit(s), taps {th}/{tl}")
        for side, gh, gl in (("p", f"bpos{k}", f"bneg{k}"), ("n", f"bneg{k}", f"bpos{k}")):
            bp = f"b{side}{k}"
            for u in range(n):
                A(f"XC{side}{k}_{u} v{side} {bp} {MIM} W={UNIT_W} L={UNIT_L}")
            A(f"XSh{side}{k} {bp} {gh} {th} vss {NFET} W={sw_w} L=0.15")
            A(f"XSl{side}{k} {bp} {gl} {tl} vss {NFET} W={sw_w} L=0.15")
            A(f"XSm{side}{k} {bp} bmid{k} vrm vss {NFET} W={sw_w} L=0.15")
    A("* comparator: StrongARM, inputs vn (+) and vp (-)")
    A(f"XMt t cmp_clk vss vss {NFET} W=4 L=0.15")
    A(f"XM1 x1 vn t vss {NFET} W=8 L=0.15 m=2")
    A(f"XM2 x2 vp t vss {NFET} W=8 L=0.15 m=2")
    A(f"XM3 op om x1 vss {NFET} W=2 L=0.15")
    A(f"XM4 om op x2 vss {NFET} W=2 L=0.15")
    A(f"XM5 op om vdd vdd {PFET} W=2 L=0.15")
    A(f"XM6 om op vdd vdd {PFET} W=2 L=0.15")
    for node, name in (("op", "7"), ("om", "8"), ("x1", "9"), ("x2", "10")):
        A(f"XMr{name} {node} cmp_clk vdd vdd {PFET} W=1 L=0.15")
    A("* output buffers: vn > vp discharges op, so cmp_p = not op")
    A(f"XIp1 cmp_p op vdd vdd {PFET} W=1 L=0.15")
    A(f"XIn1 cmp_p op vss vss {NFET} W=0.5 L=0.15")
    A(f"XIp2 cmp_n om vdd vdd {PFET} W=1 L=0.15")
    A(f"XIn2 cmp_n om vss vss {NFET} W=0.5 L=0.15")
    A(".ends")
    return "\n".join(L) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", type=pathlib.Path, default=HERE / "sarfe_ref.spice")
    args = ap.parse_args()
    args.o.write_text(netlist(), encoding="utf-8", newline="\n")
    print(f"wrote {args.o.name}: unit {C_UNIT*1e15:.2f} fF, fixed MIM "
          f"{FIX_SIDE} um square per side")


if __name__ == "__main__":
    main()
