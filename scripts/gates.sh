#!/usr/bin/env bash
# rasayan — the full gate run (port plan §6): the oracle side, the build and bundle, the parity corpus on three
# targets, benchmarks, hygiene, the lock and the audit. Adapted from ushma's scripts/gates.sh (ushma HEAD bdbcc02;
# kimiya's handoff/gates/gates.sh originally). A step whose tooling doesn't exist yet is reported as SKIP with the
# milestone that brings it, never silently passed.
#
#   bash scripts/gates.sh                                   # logs in build/gates/ (gitignored)
#   bash scripts/gates.sh > build/gates/gates.log 2>&1      # the whole run in one file
#
# Rules (CLAUDE.md, port plan §4 step 8):
# - everything is written under the repo (build/ is gitignored), never /tmp, which is a shared RAM-backed tmpfs;
# - steps run one at a time: never two wine or qemu processes at once;
# - a wine run can be truncated under load and still exit 0, so every test log must carry its summary lines, and the
#   parity summary must be the same line on every target;
# - nothing is staged: the maintainer does all git.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${OUT:-$ROOT/build/gates}"
case "$OUT" in
    /*) ;;
    *) OUT="$ROOT/$OUT" ;;
esac
case "$OUT" in
    "$ROOT"/*) ;;
    *) echo "gates.sh: OUT must be under the repo ($ROOT), not $OUT"; exit 2 ;;
esac
mkdir -p "$OUT"
cd "$ROOT" || exit 2

# Diagnostics that come from vendored lib/ (none at the pin so far). Anything else is rasayan's and fails the count.
KNOWN='^warning:lib/'

# Every step that exits nonzero is named at the end, and the run then exits 1. Each step prints its wall time.
FAILED=()
TIMES=()
T0=$(date +%s)
now() { date +%s.%N; }
elapsed() { awk -v a="$1" -v b="$2" 'BEGIN { printf "%.1f", b - a }'; }
run() { # PAT=<grep -E> [TAILN=n] run <cmd...>
    echo "=== $*"
    local s e rc
    s=$(now)
    "$@" > "$OUT/step.out" 2>&1
    rc=$?
    e=$(now)
    grep -v '^note:' "$OUT/step.out" | grep -iE "${PAT:-.}" | tail -"${TAILN:-4}"
    echo "--- exit $rc ($(elapsed "$s" "$e") s)"
    TIMES+=("$(elapsed "$s" "$e") s  $*")
    [ "$rc" -eq 0 ] || FAILED+=("$* (exit $rc)")
}
skip() { echo "=== SKIP: $*"; }

# ── 1. The Rust oracle (read-only) and the golden-value generator ────────────────────────────────────────────────
PAT='test result' TAILN=8 run cargo test --manifest-path rust-old/Cargo.toml --all-features --locked
PAT='DIFFER|STALE|missing|complete|error' run cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check
PAT='warn|error|Finished' run cargo clippy --release --locked --manifest-path scripts/oracle/Cargo.toml --all-targets -- -D warnings
PAT='FAILURES|independently checked|error' run python3 -I -B scripts/oracle/validate.py .
PAT='^ok|self-checks|STALE|DIFFER|Error|assert' run python3 -I -B scripts/oracle/mpfr_ref.py --check .
if [ -f ../ushma/src/mathx_data.cyr ]; then
    PAT='STALE|DIFFERS|current|ushma' run python3 -I -B scripts/oracle/mathx_tables.py --check-ushma
else
    PAT='STALE|DIFFERS|current' run python3 -I -B scripts/oracle/mathx_tables.py --check
fi
PAT='FAIL|checked|identical|differ' run python3 -I -B scripts/oracle/probe_tables.py . --all-targets
PAT='MISMATCH|DIFFER|cr_verify:|special values|arguments' TAILN=8 run python3 -I -B scripts/oracle/cr_verify.py --quick --targets elf,aarch64,pe

# ── 2. Build and bundle ──────────────────────────────────────────────────────────────────────────────────────────
PAT='OK|error|warn|FAIL' run cyrius build src/main.cyr build/rasayan
bw=$(grep -i warning "$OUT/step.out" | grep -viE "$KNOWN" | wc -l)
echo "rasayan warning lines (build): $bw"
[ "$bw" -eq 0 ] || FAILED+=("cyrius build ($bw warning lines)")
# [lib] only until M8 adds [lib.mcp] / [lib.ai]; --check is read-only (dist/ is committed and must be current).
PAT='.' run cyrius distlib --check

# ── 3. The parity corpus on every target, one at a time: ELF, ELF --bounds, PE (wine), aarch64 (qemu) ─────────────
for t in '' --bounds --win --aarch64; do
    log="$OUT/test${t:-_elf}.log"
    echo "=== cyrius test $t"
    s=$(now)
    cyrius test $t > "$log" 2>&1
    rc=$?
    e=$(now)
    echo "--- exit $rc ($(elapsed "$s" "$e") s)"
    TIMES+=("$(elapsed "$s" "$e") s  cyrius test $t")
    grep -E 'parity summary' "$log"
    grep -E '^[0-9]+ passed, [0-9]+ failed' "$log" | tail -1
    warns=$(grep -i warning "$log" | grep -viE "$KNOWN" | wc -l)
    echo "rasayan warning lines: $warns"
    [ "$rc" -eq 0 ] || FAILED+=("cyrius test $t (exit $rc)")
    grep -qE 'parity summary: .* 0 failing' "$log" || FAILED+=("cyrius test $t (no clean parity summary)")
    grep -qE '^[1-9][0-9]* passed, 0 failed$' "$log" || FAILED+=("cyrius test $t (no clean file summary)")
    [ "$warns" -eq 0 ] || FAILED+=("cyrius test $t ($warns rasayan warning lines)")
done
# The same corpus result on every target: one summary line, and the same per-function lines apart from the ungated
# f64_exp / f64_ln builtin shadow lines (x87 on x86_64, the polyfills on aarch64: port-parity.md §2c).
echo "=== three targets: parity summaries and per-function lines"
sums=$(for t in _elf --bounds --win --aarch64; do grep -E 'parity summary' "$OUT/test$t.log"; done | sort -u | wc -l)
echo "distinct parity summary lines across ELF, --bounds, PE, aarch64: $sums"
[ "$sums" -eq 1 ] || FAILED+=("parity summaries differ across targets ($sums distinct)")
gated() { grep -vE '^note:|~ f64_(exp|ln) builtin' "$1"; }
for t in --bounds --win --aarch64; do
    if ! diff <(gated "$OUT/test_elf.log") <(gated "$OUT/test$t.log") > "$OUT/diff$t.out"; then
        echo "test $t differs from ELF beyond the builtin shadow lines:"
        head -10 "$OUT/diff$t.out"
        FAILED+=("cyrius test $t output differs from ELF")
    fi
done

# ── 4. Benchmarks, fuzzing, coverage, hygiene, lock, audit ───────────────────────────────────────────────────────
PAT='passed|failed|ns avg' TAILN=14 run cyrius bench
grep -qE '^=== [1-9][0-9]* passed, 0 failed ===$' "$OUT/step.out" || FAILED+=("cyrius bench (no clean summary)")
PAT='passed|failed|FAIL|fuzz:' run cyrius fuzz --poison --bounds
if grep -qE '^fn ry_' src/*.cyr; then
    PAT='Functions executed|failed|below' run cyrius coverage --exec --min 80
else
    skip "cyrius coverage --exec --min 80 (no public ry_ fn in src/ until M1: nothing to measure)"
fi
echo "=== fmt/lint/doc per file"
s=$(now)
fails=0
for f in src/*.cyr tests/*.tcyr tests/*.bcyr tests/*.fcyr tests/*.cyr scripts/probes/*.cyr scripts/probes/*.bcyr docs/examples/*.cyr; do
    [ -f "$f" ] || continue
    cyrius fmt --check "$f" > /dev/null 2>&1 || { echo "fmt FAIL $f"; fails=$((fails + 1)); }
    # p12 holds the f64 traps on purpose: its diagnostics are what it measures.
    case "$f" in scripts/probes/p12_*) continue ;; esac
    cyrius lint "$f" 2>&1 | grep -E '^[1-9][0-9]* warnings?$|: warning' && { echo "  ^ lint $f"; fails=$((fails + 1)); }
done
for f in src/*.cyr docs/examples/*.cyr; do
    [ -f "$f" ] || continue
    cyrius doc --check "$f" 2>&1 | grep -E '[1-9][0-9]* undocumented' && { echo "  ^ doc $f"; fails=$((fails + 1)); }
done
e=$(now)
echo "--- fmt/lint/doc failures $fails ($(elapsed "$s" "$e") s)"
TIMES+=("$(elapsed "$s" "$e") s  fmt/lint/doc per file")
[ "$fails" -eq 0 ] || FAILED+=("fmt/lint/doc ($fails)")
PAT='.' run cyrius vet src/main.cyr
PAT='.' TAILN=2 run cyrius deps --locked
# Read by its own exit status: never piped (a pipe would report the last command's status).
PAT='ok:|FAIL|error|passed|failed' TAILN=8 run cyrius audit

# ── 5. The probes (M0 step 9): every one builds warning-free; the two compared with the oracle must agree ─────────
if [ -x scripts/probes/run.sh ]; then
    PAT='probe|FAIL|caught' TAILN=30 run bash scripts/probes/run.sh --gate
else
    skip "scripts/probes/run.sh (M0 part 3a)"
fi

# ── 6. Parity bookkeeping (from M1: every Rust public item has its Cyrius name; M9: every Rust test is ported) ────
if [ -f scripts/surface_parity.py ]; then
    PAT='MISSING|STALE|UNKNOWN|surface parity' run python3 -I -B scripts/surface_parity.py .
else
    skip "scripts/surface_parity.py (M1, with the first public fn)"
fi
if [ -f scripts/rust_test_ledger.py ]; then
    PAT='LEDGER|STALE|rust-test ledger' run python3 -I -B scripts/rust_test_ledger.py . --check
else
    skip "scripts/rust_test_ledger.py (M1)"
fi

echo "=== staged (must be empty: the maintainer does all git)"
git diff --cached --name-only
[ -z "$(git diff --cached --name-only)" ] || FAILED+=("staged files (the maintainer does all git)")

echo "=== step times"
printf '  %s\n' "${TIMES[@]}"
echo "=== total $(( $(date +%s) - T0 )) s"
if [ "${#FAILED[@]}" -gt 0 ]; then
    echo "=== FAILED: ${#FAILED[@]} step(s)"
    printf '  %s\n' "${FAILED[@]}"
    exit 1
fi
echo "=== all steps exited 0"
