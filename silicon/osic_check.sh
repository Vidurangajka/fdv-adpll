#!/usr/bin/env bash
# Sanity check of the container: the PDK the flow needs, and the tool versions.
# Run through the wrappers:  .\osic.ps1 -Check   /   ./osic.sh --check
# Exits non-zero if anything the flow depends on is missing or on the wrong PDK.
set -u
fail=0
ok()  { printf '  %-22s ok\n' "$1"; }
bad() { printf '  %-22s MISSING  (%s)\n' "$1" "$2"; fail=1; }

echo "PDK=$PDK  PDK_ROOT=$PDK_ROOT"
[ "$PDK" = sky130A ] || { echo "  PDK is $PDK, not sky130A"; fail=1; }

T="$PDK_ROOT/$PDK/libs.tech"
for f in ngspice/sky130.lib.spice magic/$PDK.magicrc netgen/${PDK}_setup.tcl; do
    [ -f "$T/$f" ] && ok "$f" || bad "$f" "$T/$f"
done

echo "versions:"
printf '  ngspice  %s\n' "$(ngspice -v 2>&1 | grep -o 'ngspice-[0-9.]*' | head -1)"
printf '  magic    %s\n' "$(magic --version 2>&1 | head -1)"
printf '  netgen   %s\n' "$(netgen -batch quit 2>&1 | grep -io 'netgen [0-9.]*' | head -1)"
printf '  klayout  %s\n' "$(klayout -v 2>&1 | head -1)"
printf '  xschem   %s\n' "$(xschem -v 2>&1 | grep -io 'xschem v[0-9.]*' | head -1)"
python3 -c "import klayout.db" 2>/dev/null && ok "klayout.db (python)" \
    || bad "klayout.db (python)" "layout generators import it"

# A one-transistor ngspice run through the real model library: proves the
# models parse, not merely that the file exists.
cat > /tmp/probe.spice <<EOF
* osic_check probe
.lib $T/ngspice/sky130.lib.spice tt
XM1 d g 0 0 sky130_fd_pr__nfet_01v8 L=1 W=5
Vd d 0 0.9
Vg g 0 0.9
.op
.control
run
let id = -i(Vd)
print id
.endc
.end
EOF
if ngspice -b /tmp/probe.spice 2>&1 | grep -q '^id = '; then ok "ngspice sky130 nfet"
else bad "ngspice sky130 nfet" "model library did not simulate"; fi

[ $fail = 0 ] && echo "PASS" || echo "FAIL"
exit $fail
