#!/usr/bin/env bash
# Process + temperature spread of the DAC's R_D * C_SAR product.
#
#   .\osic.ps1 bash ngspice/rc_spread.sh
#
# The DAC-to-ramp gain is R_D * C_SAR / T_PD, and the background calibration
# can only absorb +-3 % of it (calib.PdGainCalibration.clamp) -- so this spread
# is what the R_D trim has to cover.
#
# sky130's MOS corners (tt, ss, ff, ...) all include res_typical__cap_typical:
# resistor and capacitor process spread lives only in the separate r+c files,
# so they are loaded after the tt section here, overriding its typical
# values.  Every MOS corner run elsewhere in silicon/ has held R and C typical.
set -u
NG=$PDK_ROOT/sky130A/libs.tech/ngspice

# R: 0.35 um high-sheet poly, 2 kOhm target.  C: 1 pF MIM (m3), as C_SAR.
cat > /tmp/rc_body.spice <<'EOF'
XR a 0 0 sky130_fd_pr__res_high_po_0p35 L=2.1
Va a 0 1
XC c 0 sky130_fd_pr__cap_mim_m3_1 W=22 L=22 MF=1
Vc c 0 dc 0 ac 1
EOF

printf '%-26s %6s %10s %10s %10s\n' corner temp R_ohm C_pF RC_ns
for rc in res_typical__cap_typical res_high__cap_high res_high__cap_low \
          res_low__cap_high res_low__cap_low; do
    for t in -40 27 85; do
        cat > /tmp/rc.spice <<EOF
* rc spread
.lib "$NG/sky130.lib.spice" tt
* the tt section loaded typical R/C; the corner's .params come after and win
.include "$NG/r+c/$rc.spice"
.include "$NG/r+c/${rc}__lin.spice"
.include /tmp/rc_body.spice
.temp $t
.control
op
let r = 1 / (-i(Va))
print r
ac lin 1 1meg 1meg
let c = abs(imag(-i(Vc))) / (2 * pi * 1e6)
print c
.endc
.end
EOF
        out=$(ngspice -b /tmp/rc.spice 2>&1)
        r=$(echo "$out" | awk '/^r = /{print $3}')
        c=$(echo "$out" | awk '/^c = /{print $3}')
        if [ -z "$r" ] || [ -z "$c" ]; then
            echo "$out" | tail -15; exit 1
        fi
        awk -v n="$rc" -v t="$t" -v r="$r" -v c="$c" \
            'BEGIN{printf "%-26s %6s %10.1f %10.4f %10.4f\n", n, t, r, c*1e12, r*c*1e9}'
    done
done
