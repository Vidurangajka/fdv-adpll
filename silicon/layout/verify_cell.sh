#!/usr/bin/env bash
# Generate one cell's layout, then DRC, extract and LVS it.
#
#   .\osic.ps1 bash layout/verify_cell.sh idac_array idac_array.py idac_array_ref.spice
#
# The generic form of verify.sh (which stays as the ramp sink's record).  The
# port order is read from the reference's .subckt line, so the extracted
# netlist can be dropped into a testbench positionally.
#
# Outputs in layout/build/<cell>/ (gitignored):
#   <cell>.gds / .png / .mag, drc.txt, <cell>_lvs.spice, <cell>_pex.spice,
#   lvs.txt
#
# Exits non-zero unless DRC is clean and the circuits match uniquely.  Refuses
# an empty cell: an earlier verify.sh read the GDS wrongly and got a clean DRC
# on nothing.
set -eu
[ "${PDK:-}" = sky130A ] || { echo "PDK must be sky130A (use osic.ps1)"; exit 1; }
CELL=$1 GEN=$2 REF=$3
cd "$(dirname "$0")"
B=build/$CELL
mkdir -p "$B"
RC="$PDK_ROOT/$PDK/libs.tech/magic/$PDK.magicrc"
PORTS=$(awk -v c="$CELL" '$1==".subckt" && $2==c {for (i=3;i<=NF;i++) printf "%s ", $i}' "$REF")
[ -n "$PORTS" ] || { echo "no .subckt $CELL in $REF"; exit 1; }

python3 "$GEN" -o "$B/$CELL.gds" --png "$B/$CELL.png"

cat > "$B/run.tcl" <<EOF
gds read $CELL.gds
load $CELL
select top cell
if {[llength [cellname list allcells]] < 2 || [lindex [box values] 2] <= 0} {
    puts "EMPTY_CELL"
    quit -noprompt
}
drc style drc(full)
drc euclidean on
drc check
drc catchup
set n [drc list count total]
set fh [open drc.txt w]
puts \$fh "DRC errors: \$n"
foreach {why boxes} [drc listall why] {
    puts \$fh "\$why"
    foreach b \$boxes { puts \$fh "    \$b" }
}
close \$fh
puts "DRC_COUNT \$n"
set i 1
foreach p {$PORTS} {
    catch {port \$p make}
    if {[catch {port \$p index \$i}]} { puts "MISSING_PORT \$p" }
    incr i
}
save $CELL.mag
extract do local
extract all
ext2spice lvs
ext2spice -o ${CELL}_lvs.spice
ext2spice lvs
ext2spice cthresh 0
ext2spice -o ${CELL}_pex.spice
quit -noprompt
EOF

(cd "$B" && magic -dnull -noconsole -rcfile "$RC" run.tcl 2>&1) \
    | tee "$B/magic.log" | grep -E "DRC_COUNT|MISSING_PORT|EMPTY_CELL" || true
if grep -q EMPTY_CELL "$B/magic.log"; then echo "FAIL: empty cell"; exit 1; fi
head -30 "$B/drc.txt"

netgen -batch lvs "$B/${CELL}_lvs.spice $CELL" "$REF $CELL" \
    "$PDK_ROOT/$PDK/libs.tech/netgen/${PDK}_setup.tcl" \
    "$B/lvs.txt" > "$B/netgen.log" 2>&1 || true
grep -E "Number of devices|Number of nets|match|Mismatch|pin" "$B/lvs.txt" | tail -12

ndrc=$(head -1 "$B/drc.txt" | awk '{print $3}')
[ "$ndrc" = "0" ] || { echo "FAIL: $ndrc DRC errors"; exit 1; }
grep -q "Circuits match uniquely" "$B/lvs.txt" || { echo "FAIL: LVS"; exit 1; }
# netgen still ends on "Circuits match uniquely" when device sizes differ --
# the topology matched -- so the property check is a separate test
if grep -q "Property errors were found" "$B/lvs.txt"; then
    grep -A12 "had property errors" "$B/lvs.txt"
    echo "FAIL: LVS property errors"; exit 1
fi
echo "PASS: $CELL DRC clean, LVS match"
