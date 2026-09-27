#!/usr/bin/env bash
# Single-ended ramp (tb_ramp_se.spice) over nine PVT corners, schematic and
# extracted layout.
#
#   .\osic.ps1 bash ngspice/ramp_se.sh
#
# The extracted netlist comes from layout/verify.sh (layout/build/).
set -u
cd "$(dirname "$0")"
KEEP='sr_early|nl2 |droop_pct|v_end|casc_margin'
for net in ../layout/ramp_sink_ref.spice ../layout/build/ramp_sink_pex.spice; do
    [ -f "$net" ] || { echo "missing $net (run layout/verify.sh)"; continue; }
    echo "##### $net"
    for c in ss tt ff; do
        for t in -40 27 85; do
            sed -e "s|RAMP_NETLIST|$net|" \
                -e "s/sky130.lib.spice tt/sky130.lib.spice $c/" \
                -e "s/^\.temp 27/.temp $t/" tb_ramp_se.spice > /tmp/rse.spice
            printf '=== %s %s C: ' "$c" "$t"
            ngspice -b /tmp/rse.spice 2>&1 | grep -E "^($KEEP)" \
                | awk '{printf "%s %s  ", $1, $3}'
            echo
        done
    done
done
