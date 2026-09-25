#!/usr/bin/env bash
# Generate the ramp-sink layout, then DRC, extract and LVS it.
#
#   docker run --rm -e PDK=sky130A -v "<repo>/silicon:/work" -w /work \
#     hpretl/iic-osic-tools:latest --skip bash /work/layout/verify.sh
#
# PDK=sky130A matters: the image defaults to ihp-sg13g2.
#
# Outputs, all in layout/build/ (gitignored):
#   ramp_sink.mag         magic view of the generated GDS
#   drc.txt               magic DRC, full rule deck
#   ramp_sink_lvs.spice   extracted netlist, devices only
#   ramp_sink_pex.spice   extracted netlist with parasitic capacitance
#   lvs.txt               netgen report
#
# Exits non-zero if DRC finds anything or LVS does not match.
set -eu
export PDK=sky130A
cd "$(dirname "$0")"
mkdir -p build
RC="$PDK_ROOT/$PDK/libs.tech/magic/$PDK.magicrc"

python3 ramp_sink.py -o build/ramp_sink.gds --png build/ramp_sink.png

cat > build/run.tcl <<'EOF'
gds read ramp_sink.gds
load ramp_sink
select top cell
# A clean DRC on an empty cell is worthless, and that is exactly what a bad
# read produces -- so refuse to go on unless something was actually loaded.
if {[llength [cellname list allcells]] < 2 || [lindex [box values] 2] <= 0} {
    puts "EMPTY_CELL"
    quit -noprompt
}

# --- DRC, full deck
drc style drc(full)
drc euclidean on
drc check
drc catchup
set n [drc list count total]
set fh [open drc.txt w]
puts $fh "DRC errors: $n"
foreach {why boxes} [drc listall why] {
    puts $fh "$why"
    foreach b $boxes { puts $fh "    $b" }
}
close $fh
puts "DRC_COUNT $n"

# pin labels -> ports, in the reference's order.  The GDS carries both a text
# label and a pin shape for each, so address them by name, not by cursor box.
set i 1
foreach p {vdd vss pre outp nbias cbias} {
    # the pin shapes usually arrive as ports already; `make` covers the case
    # where they do not, and `index` fixes the order either way -- the PEX
    # testbench instantiates the subcircuit positionally.
    catch {port $p make}
    if {[catch {port $p index $i}]} { puts "MISSING_PORT $p" }
    incr i
}
puts "PORTS=[port last]"
save ramp_sink.mag

# --- extraction: devices only for LVS, then with capacitance for simulation
extract do local
extract all
ext2spice lvs
ext2spice -o ramp_sink_lvs.spice
ext2spice lvs
ext2spice cthresh 0
ext2spice -o ramp_sink_pex.spice
quit -noprompt
EOF

(cd build && magic -dnull -noconsole -rcfile "$RC" run.tcl 2>&1) \
    | tee build/magic.log | grep -E "DRC_COUNT|MISSING_PORT|PORTS=|EMPTY_CELL|rror" || true
if grep -q EMPTY_CELL build/magic.log; then echo "FAIL: GDS read produced an empty cell"; exit 1; fi

head -40 build/drc.txt

netgen -batch lvs \
    "build/ramp_sink_lvs.spice ramp_sink" \
    "ramp_sink_ref.spice ramp_sink" \
    "$PDK_ROOT/$PDK/libs.tech/netgen/${PDK}_setup.tcl" \
    build/lvs.txt > build/netgen.log 2>&1 || true
tail -n 30 build/lvs.txt

ndrc=$(head -1 build/drc.txt | awk '{print $3}')
if [ "$ndrc" != "0" ]; then echo "FAIL: $ndrc DRC errors"; exit 1; fi
if ! grep -q "Circuits match uniquely" build/lvs.txt; then echo "FAIL: LVS"; exit 1; fi
echo "PASS: DRC clean, LVS match"
