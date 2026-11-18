#!/usr/bin/env bash
# rasayan — the M0 step 9 probes (port plan §4 step 9; write-ups in docs/development/port-parity.md §3).
#
#   bash scripts/probes/run.sh            # build every probe, run each on x86_64 ELF and print what it measures
#   bash scripts/probes/run.sh --gate     # what scripts/gates.sh runs (below)
#
# Both modes first write the references in Python (scripts/probes/probe_substeps.py, probe_solver.py, probe_mcp.py,
# probe_json.py: rust-old/src read as text, nothing run in Rust), into build/probe-substeps/, build/probe-solver/,
# build/probe-mcp/ and build/probe-json/.
# Then:
# - every probe builds on x86_64 with no warning line, except p12, which holds the f64 traps on purpose and must give
#   exactly its 8 diagnostics; p12's exponent-literal trap (-D P12_EXPONENT) must be a compile error;
# - p1 (sub-step counts) and p4 (the implicit NADH solves) must agree with the reference on every row, on x86_64 ELF,
#   aarch64 (qemu) and x86_64 PE (wine), one target at a time, each log checked for its summary line (wine can truncate
#   output under load and still exit 0); p4's control build (-D P4_ALT, the algebraically equal secant step) must NOT
#   agree, which shows the comparison sees a 1-ulp difference;
# - p2 (copies without the heap), p3 (wide results) and p5 (the network's 10^6-step soak, M5) check their own
#   claims and exit 1 if one fails;
# - p6 (bote 3.8.15, part 3b) answers the reference MCP requests on x86_64 ELF, aarch64 and PE, and p6_compare.py must
#   find every tools/call payload and the inputSchema byte-identical and no difference but the recorded ones;
# - p7 (hoosh 2.12.3, P8; since M8 part 2 rasayan's ai client itself): the reference client probe_ai.py and p7_hoosh
#   each make the client's eight calls against their own
#   loopback stub (p7_hoosh_stub.py, ports 18741 and 18742), and p7_compare.py must find the request lines, bodies,
#   Content-Type and Accept, and every outcome line identical (x86_64 ELF);
# - p10 (the manifest and the bundles) runs p10_manifest.sh, which checks its own claims (consumers on three targets);
# - p11 (JSON, P7) must reproduce every serde_json document byte for byte and read every one bayan accepts back
#   exactly, and fmtf must match serde_json and Display on every double (x86_64 ELF and PE; aarch64 outside --gate);
#   its control build (-D P11_ALT, bayan's own printer) must NOT agree;
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

# The references (Python; rust-old/src read as text).
for ref in probe_substeps probe_solver probe_mcp probe_json; do
    if python3 -I -B "scripts/probes/$ref.py" > "$OUT/$ref.log" 2>&1; then
        grep -E "^$ref:" "$OUT/$ref.log"
    else
        cat "$OUT/$ref.log"
        FAILED+=("reference $ref")
    fi
done

# Every probe builds; only p12 warns, and exactly as often as it has diagnosed traps.
for p in p0_math p1_substeps p2_memory p3_wide_results p4_solver p5_soak p6_bote p7_hoosh p8_span_cost p11_json \
    p12_f64_hygiene; do
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
# A wine run that exits 0 without its summary line was truncated (wine drops buffered output under load: port-parity
# §2f); it is run once more, and says so. Nothing else is retried, and the line must then be there.
runs() {
    local label="$1" want="$2"
    shift 2
    "$@" > "$OUT/$label.log" 2>&1
    local rc=$?
    if [ "$rc" -eq 0 ] && [ "$1" = wine ] && ! grep -qE "$want" "$OUT/$label.log"; then
        echo "$label: wine output truncated ($(wc -l < "$OUT/$label.log") lines); run once more"
        "$@" > "$OUT/$label.log" 2>&1
        rc=$?
    fi
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
# p5, the network's 10^6-step soak (M5 part 3), checks its own claims.
runs p5_soak 'probe 5: .*: ok$' "./$OUT/p5_soak"

# p6 on the three targets: bote's answers to the oracle's requests, compared by p6_compare.py.
p6cmp() { # <label> <cmd...>: runs p6, keeps its answers, compares them with the oracle's
    local label="$1"
    shift
    "$@" 2> "$OUT/$label.err" | tr -d '\r' > "$OUT/$label.tsv"
    local rc=${PIPESTATUS[0]}
    [ "$rc" -eq 0 ] || FAILED+=("$label (exit $rc)")
    runs "$label" 'probe 6: unexpected differences 0$' python3 -I -B scripts/probes/p6_compare.py \
        build/probe-mcp/rust.tsv "$OUT/$label.tsv"
}
p6cmp p6_bote-elf "./$OUT/p6_bote"
if cyrius build --aarch64 scripts/probes/p6_bote.cyr "$OUT/p6_bote-aarch64" > "$OUT/p6_bote-aarch64.build" 2>&1; then
    p6cmp p6_bote-aarch64 qemu-aarch64 "./$OUT/p6_bote-aarch64"
else
    FAILED+=("build --aarch64 p6_bote")
fi
if cyrius build --win scripts/probes/p6_bote.cyr "$OUT/p6_bote.exe" > "$OUT/p6_bote-pe.build" 2>&1; then
    p6cmp p6_bote-pe wine "./$OUT/p6_bote.exe"
else
    FAILED+=("build --win p6_bote")
fi

# p7: the reference client and the Cyrius requests, each against its own hoosh 2.12.3 stub (8 requests each).
rm -f "$OUT/p7-ref.log" "$OUT/p7-cyrius.log"
python3 -I -B scripts/probes/p7_hoosh_stub.py 18741 "$OUT/p7-ref.log" 8 &
stub_r=$!
python3 -I -B scripts/probes/p7_hoosh_stub.py 18742 "$OUT/p7-cyrius.log" 8 &
stub_c=$!
python3 -I -B -c 'import time; time.sleep(0.5)'
python3 -I -B scripts/probes/probe_ai.py http://127.0.0.1:18741 > "$OUT/p7-ref.out" 2>&1 || FAILED+=("reference probe_ai")
"./$OUT/p7_hoosh" http://127.0.0.1:18742 > "$OUT/p7-cyrius.out" 2>&1 || FAILED+=("p7_hoosh")
wait "$stub_r" "$stub_c" 2> /dev/null
runs p7_hoosh 'probe 7: differences 0$' python3 -I -B scripts/probes/p7_compare.py "$OUT/p7-ref.log" \
    "$OUT/p7-cyrius.log" "$OUT/p7-ref.out" "$OUT/p7-cyrius.out"

# p10: the manifest, the bundles and their consumers (its own checks).
runs p10_manifest 'probe 10: all checks passed$' bash scripts/probes/p10_manifest.sh

# p11 on x86_64 ELF and PE (fmtf is integer arithmetic; aarch64 outside --gate), and its control.
p11t="elf pe"
[ "$GATE" -eq 0 ] && p11t="elf pe aarch64"
for t in $p11t; do
    case $t in
        elf) runs p11_json-elf 'probe 11: failed checks 0$' "./$OUT/p11_json" ;;
        pe)
            if cyrius build --win scripts/probes/p11_json.cyr "$OUT/p11_json.exe" > "$OUT/p11_json-pe.build" 2>&1; then
                runs p11_json-pe 'probe 11: failed checks 0$' wine "./$OUT/p11_json.exe"
            else
                FAILED+=("build --win p11_json")
            fi
            ;;
        aarch64)
            if cyrius build --aarch64 scripts/probes/p11_json.cyr "$OUT/p11_json-aarch64" > "$OUT/p11_json-aarch64.build" 2>&1; then
                runs p11_json-aarch64 'probe 11: failed checks 0$' qemu-aarch64 "./$OUT/p11_json-aarch64"
            else
                FAILED+=("build --aarch64 p11_json")
            fi
            ;;
    esac
done
if cyrius build -D P11_ALT scripts/probes/p11_json.cyr "$OUT/p11alt" > "$OUT/p11alt.build" 2>&1; then
    if "./$OUT/p11alt" > "$OUT/p11alt.log" 2>&1; then
        FAILED+=("p11 control (-D P11_ALT) was not caught")
    else
        echo "p11 control (-D P11_ALT): caught, as it must be: $(grep -E 'probe 11: writer differs' "$OUT/p11alt.log")"
    fi
else
    FAILED+=("build -D P11_ALT p11_json")
fi

if [ "$GATE" -eq 0 ]; then
    for p in p1_substeps-elf p4_solver-elf p2_memory p3_wide_results p6_bote-elf p7_hoosh p11_json-elf; do
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
