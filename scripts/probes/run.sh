#!/usr/bin/env bash
# rasayan — the M0 step 9 probes (port plan §4 step 9; write-ups in docs/development/port-parity.md §3).
#
#   bash scripts/probes/run.sh            # build every probe, run each on x86_64 ELF and print what it measures
#   bash scripts/probes/run.sh --gate     # what scripts/gates.sh runs (below)
#
# Both modes first run the oracle's probe bins (scripts/oracle/src/bin/probe_substeps.rs, probe_solver/), which write
# build/probe-substeps/ and build/probe-solver/. Then:
# - every probe builds on x86_64 with no warning line, except p12, which holds the f64 traps on purpose and must give
#   exactly its 8 diagnostics; p12's exponent-literal trap (-D P12_EXPONENT) must be a compile error;
# - p1 (sub-step counts) and p4 (the implicit NADH solves) must agree with the oracle on every row, on x86_64 ELF,
#   aarch64 (qemu) and x86_64 PE (wine), one target at a time, each log checked for its summary line (wine can truncate
#   output under load and still exit 0); p4's control build (-D P4_ALT, the algebraically equal secant step) must NOT
#   agree, which shows the comparison sees a 1-ulp difference;
# - p2 (copies without the heap) and p3 (wide results) check their own claims and exit 1 if one fails;
# - --gate leaves out the pure measurements (p0, p8, p12's run), which only build there.
# Everything is written under build/ (gitignored), never /tmp.
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT" || exit 2
OUT=build/probes
mkdir -p "$OUT"
GATE=0
[ "${1:-}" = "--gate" ] && GATE=1
FAILED=()

# The oracle's side.
for bin in probe_substeps probe_solver; do
    if cargo run -q --release --locked --manifest-path scripts/oracle/Cargo.toml --bin "$bin" > "$OUT/$bin.log" 2>&1; then
        grep -E "^$bin:" "$OUT/$bin.log"
    else
        cat "$OUT/$bin.log"
        FAILED+=("oracle bin $bin")
    fi
done

# Every probe builds; only p12 warns, and exactly as often as it has diagnosed traps.
for p in p0_math p1_substeps p2_memory p3_wide_results p4_solver p8_span_cost p12_f64_hygiene; do
    if ! cyrius build "scripts/probes/$p.cyr" "$OUT/$p" > "$OUT/$p.build" 2>&1; then
        FAILED+=("build $p")
        grep -vE '^note' "$OUT/$p.build" | head -5
        continue
    fi
    w=$(grep -c 'warning' "$OUT/$p.build")
    expect=0
    [ "$p" = p12_f64_hygiene ] && expect=8
    echo "probe build $p: $w warning lines (expected $expect)"
    [ "$w" -eq "$expect" ] || FAILED+=("build $p: $w warning lines, expected $expect")
done
if cyrius build -D P12_EXPONENT scripts/probes/p12_f64_hygiene.cyr "$OUT/p12e" > "$OUT/p12e.build" 2>&1; then
    FAILED+=("p12: an exponent literal compiled")
else
    grep -qE "expected ';', got identifier 'e'" "$OUT/p12e.build" || FAILED+=("p12: the exponent literal's error changed")
    echo "probe build p12 -D P12_EXPONENT: a compile error, as expected"
fi

# runs <label> <summary regex> <cmd...>: runs one probe binary, keeps its log, requires exit 0 and the summary line.
runs() {
    local label="$1" want="$2"
    shift 2
    "$@" > "$OUT/$label.log" 2>&1
    local rc=$?
    printf '%s: %s\n' "$label" "$(grep -E "$want" "$OUT/$label.log" | tail -1)"
    [ "$rc" -eq 0 ] || FAILED+=("$label (exit $rc)")
    grep -qE "$want" "$OUT/$label.log" || FAILED+=("$label (no summary line)")
}

# p1 and p4 on the three targets, one at a time.
for p in p1_substeps p4_solver; do
    want='probe 1: idiom mismatches 0$'
    [ "$p" = p4_solver ] && want='probe 4: mismatching rows 0$'
    runs "$p-elf" "$want" "./$OUT/$p"
    if cyrius build --aarch64 "scripts/probes/$p.cyr" "$OUT/$p-aarch64" > "$OUT/$p-aarch64.build" 2>&1; then
        runs "$p-aarch64" "$want" qemu-aarch64 "./$OUT/$p-aarch64"
    else
        FAILED+=("build --aarch64 $p")
    fi
    if cyrius build --win "scripts/probes/$p.cyr" "$OUT/$p.exe" > "$OUT/$p-pe.build" 2>&1; then
        runs "$p-pe" "$want" wine "./$OUT/$p.exe"
    else
        FAILED+=("build --win $p")
    fi
done
# p4's control: the textbook secant step must be caught.
if cyrius build -D P4_ALT scripts/probes/p4_solver.cyr "$OUT/p4alt" > "$OUT/p4alt.build" 2>&1; then
    "./$OUT/p4alt" > "$OUT/p4alt.log" 2>&1
    rc=$?
    if [ "$rc" -eq 0 ]; then
        FAILED+=("p4 control (-D P4_ALT) was not caught")
    else
        echo "p4 control (-D P4_ALT): caught, as it must be: $(grep -E 'probe 4: mismatching rows' "$OUT/p4alt.log")"
    fi
else
    FAILED+=("build -D P4_ALT p4_solver")
fi

# p2 and p3 check their own claims.
runs p2_memory 'probe 2: .*: ok$' "./$OUT/p2_memory"
runs p3_wide_results 'probe 3: .*: ok$' "./$OUT/p3_wide_results"

if [ "$GATE" -eq 0 ]; then
    for p in p1_substeps-elf p4_solver-elf p2_memory p3_wide_results; do
        echo "--- $p"
        grep -v 'per op in ps' "$OUT/$p.log"
    done
    for p in p8_span_cost p12_f64_hygiene p0_math; do
        echo "--- $p"
        [ "$p" = p12_f64_hygiene ] && grep -E 'warning' "$OUT/$p.build"
        "./$OUT/$p" 2>&1 | grep -v 'per op in ps'
    done
fi

if [ "${#FAILED[@]}" -gt 0 ]; then
    echo "probes: FAILED ${#FAILED[@]}"
    printf '  %s\n' "${FAILED[@]}"
    exit 1
fi
echo "probes: all checks passed"
