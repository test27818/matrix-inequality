#!/bin/sh
# Rebuild and replay one documented family of mathematical certificates.
# Does not claim uncovered domains, or substitute for the analytical proofs.
set -eu
cd "$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
export PYTHONDONTWRITEBYTECODE=1
tmp=$(mktemp -d "${TMPDIR:-/tmp}/psd-math-replay.XXXXXXXX")
trap 'rm -rf "$tmp"' EXIT HUP INT TERM
CXX=${CXX:-c++}
compile() {
    "$CXX" -std=c++20 -O2 -ffp-contract=off -fno-fast-math -fsanitize=undefined         -fno-sanitize-recover=all "$@"
}
py() { script=$1; shift; python3 "certificates/$script" "$@"; }
shallow() {
    py run_signed_curvature.py 4 4 8 16 6
    py run_signed_curvature.py 6 3 8 16 6
    compile certificates/short_zero_residual.cpp -o "$tmp/short-zero"
    "$tmp/short-zero" 0 3 6 8 64 22 80000
    "$tmp/short-zero" 0 3 8 24 32 19 40000
    py zero_cross_helmholtz.py --root all
    py verify_zero_far48.py
}
mixed() {
    py compact_cross_helmholtz.py 3 9 24 32 3/4 6
    py compact_cross_helmholtz.py 9 12 24 32 7/8 6
    py compact_cross_helmholtz.py 12 16 24 32 29/32 7
    py shallow_cross_helmholtz.py --root all --mesh 6
    py far_shallow_bridge.py --mesh 6
    py verify_far_shallow_bridge_constants.py
    py verify_far_shallow_physical.py
    # That experimental checker does not signal an uncertified root by
    # exit status; require the proof status field for each exact root.
    for j in 0 1; do
        py compact_separated_helmholtz.py --root "$j" >"$tmp/compact-$j"
        cat "$tmp/compact-$j"
        awk -v root="$j" '$1==root && $6=="CERTIFIED"{count++} END{exit count!=1}' "$tmp/compact-$j"
    done
}
moderate() {
    py moderate_helmholtz.py --root all --mesh 6
    py outer_moderate_helmholtz.py --mesh 6
    py verify_outer_moderate_constants.py
    py compact_gap_helmholtz.py --mesh 6
    py compact_near_errors.py
    compile certificates/compact_near_spectrum.cpp -Wl,-l:libmpfr.so.6 -o "$tmp/near"
    "$tmp/near" 16 24 0 4 128 18 4
    py compact_mid_errors.py
    compile certificates/compact_mid_spectrum.cpp -Wl,-l:libmpfr.so.6 -o "$tmp/mid12"
    "$tmp/mid12" 12 16 0 6 256 18 4
    py compact_continuum_mixed_bounds.py
    compile certificates/compact_continuum_mixed.cpp -o "$tmp/mid8"
    "$tmp/mid8" 8 12 8 16 2 64 5 8
    "$tmp/mid8" 3 8 3 16 2 64 7 16
    compile certificates/compact_angle_residual.cpp -o "$tmp/angle"
    "$tmp/angle" 8 12 16 24 32 24 80000
}
deep() {
    py verify_deep_separation.py
    sh certificates/replay_deep_continuum.sh
}
case "${1:-all}" in
    shallow) shallow ;;
    mixed) mixed ;;
    moderate) moderate ;;
    deep) deep ;;
    all) py verify_cover_geometry.py; shallow; mixed; moderate; deep ;;
    *) echo "Usage: $0 [all|shallow|mixed|moderate|deep]" >&2; exit 2 ;;
esac
