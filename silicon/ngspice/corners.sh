#!/usr/bin/env bash
# Sweep a ramp testbench over sky130 process corners and temperature.
#
#   docker run --rm -v "<repo>/silicon:/work" -w /work \
#     hpretl/iic-osic-tools:latest --skip bash /work/ngspice/corners.sh
#
# Takes the testbench path as $1 (default ngspice/tb_ramp_casc.spice).
#
# The testbench pins its corner with a literal `.lib ... tt`, so the corner is
# swapped by rewriting that line into a scratch copy rather than by parameter --
# ngspice has no way to select a .lib section from the command line.
set -u

TB="${1:-ngspice/tb_ramp_casc.spice}"

# Only the measurements that decide whether the cascode holds.  .meas PARAM
# results come back inside ngspice's "insertnumber: fails" diagnostic rather
# than on their own line -- that message is cosmetic (the value is computed and
# printed, only the back-substitution into the netlist copy fails), so the
# numbers are recovered from it here.
KEEP='sr_early_c param|droop_pct_c param|nl2_c param|bias_headroom param|mirror_match param|casc_margin param'

for c in ss tt ff; do
    for t in -40 27 85; do
        sed -e "s/sky130.lib.spice tt/sky130.lib.spice ${c}/" \
            -e "s/^\.tran 2p 12n/.temp ${t}\n.tran 2p 12n/" \
            "$TB" > /tmp/corner.spice
        printf '=== %s @ %s C ===\n' "$c" "$t"
        ngspice -b /tmp/corner.spice 2>&1 \
            | grep -E "$KEEP" \
            | sed -e 's/.*meas tran //' -e 's/ param=/ = /' -e 's/".*//'
    done
done
