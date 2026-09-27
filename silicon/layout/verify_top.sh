#!/usr/bin/env bash
# The whole detector: generate, DRC, extract, LVS.
#
#   .\osic.ps1 bash layout/verify_top.sh
#
# Needs the four analog cells verified first (verify_cell.sh writes their GDS
# into build/<cell>/) and a finished LibreLane run of digital/.  The reference
# is fdvpd_ref.v (top_netlist.py) over the cells' *_ref.spice, the macro's
# powered netlist and the hd library's SPICE.
#
# Outputs in layout/build/fdvpd/: fdvpd.gds / .png, drc.txt, fdvpd_lvs.spice,
# lvs.txt.  Exits non-zero unless DRC is clean and LVS matches uniquely.
set -eu
[ "${PDK:-}" = sky130A ] || { echo "PDK must be sky130A (use osic.ps1)"; exit 1; }
cd "$(dirname "$0")"
B=build/fdvpd
mkdir -p "$B"
# the newest run that finished (top.py takes the same one)
RUN=$(for r in $(ls -d ../digital/runs/RUN_* | sort -r); do
          [ -f "$r/final/pnl/fdvpd_ctrl.pnl.v" ] && { echo "$r"; break; }; done)
[ -n "$RUN" ] || { echo "no finished LibreLane run in ../digital/runs"; exit 1; }
echo "macro: $RUN"
RC="$PDK_ROOT/$PDK/libs.tech/magic/$PDK.magicrc"
HD="$PDK_ROOT/$PDK/libs.ref/sky130_fd_sc_hd/spice/sky130_fd_sc_hd.spice"

python3 top.py -o "$B/fdvpd.gds" --png "$B/fdvpd.png"
python3 top_netlist.py
python3 bias.py --ref build/bias/bias_ref_sized.spice

cat > "$B/run.tcl" <<EOF
gds read fdvpd.gds
load fdvpd
select top cell
drc style drc(full)
drc euclidean on
drc check
drc catchup
set n [drc list count total]
set fh [open drc.txt w]
puts \$fh "DRC errors: \$n"
foreach {why boxes} [drc listall why] {
    puts \$fh "\$why"
    foreach b [lrange \$boxes 0 9] { puts \$fh "    \$b" }
}
close \$fh
puts "DRC_COUNT \$n"
extract do local
extract all
ext2spice lvs
ext2spice -o fdvpd_lvs.spice
quit -noprompt
EOF
(cd "$B" && magic -dnull -noconsole -rcfile "$RC" run.tcl 2>&1) \
    | tee "$B/magic.log" | grep -E "DRC_COUNT" || true
head -40 "$B/drc.txt"

cat > "$B/lvs.tcl" <<EOF
set l [readnet spice $B/fdvpd_lvs.spice]
set r [readnet spice $HD]
readnet spice idac_ref.spice \$r
readnet spice ramp_sink_ref.spice \$r
readnet spice sarfe_ref.spice \$r
readnet spice build/bias/bias_ref_sized.spice \$r
readnet verilog $RUN/final/pnl/fdvpd_ctrl.pnl.v \$r
readnet verilog fdvpd_ref.v \$r
lvs "\$l fdvpd" "\$r fdvpd" $PDK_ROOT/$PDK/libs.tech/netgen/${PDK}_setup.tcl $B/lvs.txt
EOF
netgen -batch source "$B/lvs.tcl" > "$B/netgen.log" 2>&1 || true
grep -E "Circuits match|do not match|Property errors|Netlists do not" "$B/lvs.txt" | tail -8

ndrc=$(head -1 "$B/drc.txt" | awk '{print $3}')
[ "$ndrc" = "0" ] || { echo "FAIL: $ndrc DRC errors"; exit 1; }
tail -3 "$B/lvs.txt" | grep -q "Circuits match uniquely" || { echo "FAIL: LVS"; exit 1; }
if grep -q "Property errors were found" "$B/lvs.txt"; then
    echo "FAIL: LVS property errors"; exit 1
fi
echo "PASS: fdvpd DRC clean, LVS match"
