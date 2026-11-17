#!/usr/bin/env bash
# Probe 10 (port plan P5, §3.2, §4 step 9; port-parity.md §3): the manifest and the three bundles, measured in a
# scratch copy (build/p10/, gitignored) so nothing tracked changes. Adapted from ushma's scripts/probes/p10_manifest.sh
# (ushma commit 34e53b2).
#
# 1. The tracked manifest as it is: sakshi and log vendored and locked, bote (optional) locked with its graph,
#    `cyrius deps --locked` clean, and what [lib], [lib.mcp] and [lib.ai] bundle today.
# 2. The scratch copy adds a stand-in to [lib.mcp] (p10_manifest/mcp.cyr: one tool on bote-core, probe 6's shape)
#    and to [lib.ai] (p10_manifest/ai.cyr: a GET on sandhi, probe 7's shape); `cyrius distlib --all`, the .deps
#    sidecars and requires blocks, `distlib --check`, and whether dist/rasayan.cyr names anything of bote or sandhi.
# 3. Consumers of each bundle, as a downstream project declares them (a path dep on the scratch copy, the sidecar's
#    leaves as stdlib; the mcp one also declares bote-core itself): what `cyrius deps` vendors for each (a core
#    consumer must get neither bote nor sandhi), then a build for x86_64 ELF, PE and aarch64 and a run (wine, qemu;
#    one at a time), the ai one against the probe-7 stub on 127.0.0.1:18743.
#
#   bash scripts/probes/p10_manifest.sh            # prints the record; exits 1 if a check fails
set -u
here=$(cd "$(dirname "$0")/../.." && pwd)
OUT=$here/build/p10
rm -rf "$OUT"
mkdir -p "$OUT"
cd "$here" || exit 2
FAILED=()
# Diagnostics that are not rasayan's: from vendored lib/ (PE syscall routing in dynlib/fdlopen, filed by ushma) and
# sigil's TLS buffers in sandhi's stack (large static data, prakash's record).
KNOWN='warning:lib/|large static data'
say() { printf '\n=== %s\n' "$*"; }
quiet() { grep -v -E '^\s*lib/|rm <file>|^note:|outside the declared footprint|^\s*$'; }

say "1. the tracked manifest"
for leaf in sakshi log; do
    if [ -f "lib/$leaf.cyr" ] && grep -q " lib/$leaf.cyr\$" cyrius.lock; then
        echo "lib/$leaf.cyr: vendored and locked ($(grep " lib/$leaf.cyr\$" cyrius.lock | cut -c1-16)…)"
    else
        echo "lib/$leaf.cyr: NOT vendored and locked"
        FAILED+=("$leaf not vendored and locked")
    fi
done
grep '^commit' cyrius.lock | awk '{printf "locked: %s %s @%s\n", $3, $5, substr($2, 1, 7)}'
grep -q $'^commit\t[0-9a-f]*\tbote\t' cyrius.lock || FAILED+=("bote not locked")
if cyrius deps --locked > "$OUT/locked.log" 2>&1; then
    echo "cyrius deps --locked: $(grep -- '--locked' "$OUT/locked.log" | tail -1)"
else
    tail -3 "$OUT/locked.log"
    FAILED+=("cyrius deps --locked")
fi
for b in dist/rasayan.cyr dist/rasayan-mcp.cyr dist/rasayan-ai.cyr; do
    printf '%s: %s bytes, modules: %s\n' "$b" "$(wc -c < "$b")" "$(grep -c '^# --- ' "$b")"
done

say "2. the scratch copy: [lib.mcp] and [lib.ai] with their stand-ins"
S=$OUT/rasayan
mkdir -p "$S"
cp -r cyrius.cyml cyrius.lock VERSION src lib dist "$S/"
cp scripts/probes/p10_manifest/mcp.cyr "$S/src/p10_mcp.cyr"
cp scripts/probes/p10_manifest/ai.cyr "$S/src/p10_ai.cyr"
python3 -I -B - "$S/cyrius.cyml" <<'EOF'
import sys
p = sys.argv[1]
s = open(p).read()
mcp_at = s.index("\n[lib.mcp]\n")
ai_at = s.index("\n[lib.ai]\n")
mcp = s[mcp_at:ai_at].replace('    "src/membrane_data.cyr",\n]', '    "src/membrane_data.cyr",\n    "src/p10_mcp.cyr",\n]')
ai_end = s.index("\n[test]\n")
ai = s[ai_at:ai_end].replace('    "src/ai_data.cyr",\n]', '    "src/ai_data.cyr",\n    "src/p10_ai.cyr",\n]')
s = s[:mcp_at] + mcp + ai + s[ai_end:]
assert "p10_mcp" in s and "p10_ai" in s
open(p, "w").write(s)
EOF
(cd "$S" && cyrius distlib --all 2>&1 | quiet | sed "s|$S/||g")
for d in "$S"/dist/*.deps; do
    printf '%s: %s\n' "$(basename "$d")" "$(grep -v '^#' "$d" | tr '\n' ' ')"
done
for b in rasayan rasayan-mcp rasayan-ai; do
    printf '%s.cyr: %s bytes; requires: %s\n' "$b" "$(wc -c < "$S/dist/$b.cyr")" \
        "$(sed -n '/^# Requires/,/^$/p' "$S/dist/$b.cyr" | grep '^include' | sed 's/include "lib\/\(.*\)\.cyr"/\1/' | tr '\n' ' ')"
done
names=$(grep -c -E 'dispatcher_|codec_process|tool_registry|sandhi_' "$S/dist/rasayan.cyr")
echo "dist/rasayan.cyr lines naming bote or sandhi: $names"
[ "$names" -eq 0 ] || FAILED+=("the core bundle names bote or sandhi")
if (cd "$S" && cyrius distlib --check > "$OUT/check.log" 2>&1); then
    echo "cyrius distlib --check: $(grep -c 'current:' "$OUT/check.log") bundles current"
else
    quiet < "$OUT/check.log" | tail -5
    FAILED+=("distlib --check in the scratch copy")
fi
echo "lib/sandhi.cyr in the scratch copy: $([ -f "$S/lib/sandhi.cyr" ] && echo present || echo 'absent (the [lib.ai] verify compiled the toolchain'"'"'s)')"

# consumer <name> <bundle> <extra manifest text> <define or empty>
consumer() {
    local name=$1 bundle=$2 extra=$3 def=$4 c=$OUT/$1
    mkdir -p "$c/src"
    cp scripts/probes/p10_manifest/consumer.cyr "$c/src/main.cyr"
    local side leaves
    side=$S/dist/${bundle%.cyr}.deps
    leaves=$(grep -v '^#' "$side" 2>/dev/null | sed 's/.*/"&"/' | paste -sd, -)
    cat > "$c/cyrius.cyml" <<EOF
[package]
name = "$name"
version = "0.0.1"
language = "cyrius"
cyrius = "6.7.18"

[build]
entry = "src/main.cyr"
output = "build/$name"
dce = true

[deps]
stdlib = ["string", "fmt", "alloc", "syscalls", "io"${leaves:+, $leaves}]

[deps.rasayan]
path = "../rasayan"
modules = ["dist/$bundle"]
$extra
EOF
    say "3. consumer $name (dist/$bundle; stdlib leaves from its sidecar: ${leaves:-none})"
    (cd "$c" && cyrius deps > "$c/deps.log" 2>&1) || FAILED+=("$name: cyrius deps")
    printf 'vendored: %s\n' "$(cd "$c/lib" && ls | grep -E '^(rasayan|bote|libro|majra|patra|sandhi|tls|sigil|sakshi|log)' | tr '\n' ' ')"
    printf 'locked commits: %s\n' "$(grep '^commit' "$c/cyrius.lock" 2>/dev/null | awk '{print $3}' | tr '\n' ' ')"
    local t
    for t in elf win aarch64; do
        local flag="" ext=""
        [ "$t" = win ] && flag=--win && ext=.exe
        [ "$t" = aarch64 ] && flag=--aarch64
        local log rc
        # shellcheck disable=SC2086
        log=$(cd "$c" && cyrius build $flag $def src/main.cyr "build/$name-$t$ext" 2>&1)
        rc=$?
        local warns
        warns=$(printf '%s\n' "$log" | grep -i 'warning' | grep -v -E "$KNOWN" | wc -l)
        printf '%-8s build rc %d, %s bytes, %d warning lines (known, from vendored lib/ and TLS: %d)\n' "$t" "$rc" \
            "$(wc -c < "$c/build/$name-$t$ext" 2>/dev/null || echo 0)" "$warns" \
            "$(printf '%s\n' "$log" | grep -i 'warning' | grep -c -E "$KNOWN")"
        if [ "$name-$t" = c_mcp-win ]; then
            # Upstream (port-parity.md §3, probe 10; docs/development/upstream/bote/): a bote-core consumer compiles
            # bote's whole graph, and lib/majra.cyr has no PE futex. Expected to fail until that is fixed.
            if [ "$rc" -ne 0 ] && printf '%s\n' "$log" | grep -q "error:lib/majra.cyr:.*undefined variable 'SYS_FUTEX'"; then
                echo "win      expected failure: lib/majra.cyr (bote's graph) has no SYS_FUTEX on PE (recorded constraint of the -mcp bundle)"
                continue
            fi
            [ "$rc" -eq 0 ] && echo "win      builds now: the bote/majra PE note can be retired"
        fi
        [ "$rc" -eq 0 ] || { FAILED+=("$name: build $t"); continue; }
        [ "$warns" -eq 0 ] || FAILED+=("$name: build $t warns")
        local run rrc
        case $t in
            elf) run=$("$c/build/$name-$t" 2>/dev/null) ;;
            win) run=$(wine "$c/build/$name-$t$ext" 2>/dev/null | tr -d '\r') ;;
            aarch64) run=$(qemu-aarch64 "$c/build/$name-$t" 2>/dev/null) ;;
        esac
        rrc=$?
        printf '%-8s run rc %d: %s\n' "$t" "$rrc" "$(printf '%s' "$run" | tr '\n' ' ')"
        printf '%s' "$run" | grep -q '^core 4613303445314885481$' || FAILED+=("$name $t: core answer")
        case $def in
            *P10_MCP*) printf '%s' "$run" | grep -q '^mcp 1$' || FAILED+=("$name $t: mcp answer") ;;
            *P10_AI*) printf '%s' "$run" | grep -q '^ai 200$' || FAILED+=("$name $t: ai answer") ;;
        esac
    done
}

consumer c_core rasayan.cyr "" ""
if ls "$OUT/c_core/lib" | grep -q -E '^(bote|sandhi|libro|majra)'; then FAILED+=("c_core vendored bote or sandhi"); fi
consumer c_mcp rasayan-mcp.cyr '
[deps.bote]
git = "https://github.com/MacCracken/bote.git"
tag = "3.8.15"
modules = ["dist/bote-core.cyr"]' "-D P10_MCP"
python3 -I -B scripts/probes/p7_hoosh_stub.py 18743 "$OUT/stub.log" 3 &
stub=$!
python3 -I -B -c 'import time; time.sleep(0.5)'
consumer c_ai rasayan-ai.cyr "" "-D P10_AI"
kill $stub 2> /dev/null
wait $stub 2> /dev/null

if [ "${#FAILED[@]}" -gt 0 ]; then
    say "probe 10: FAILED ${#FAILED[@]}"
    printf '  %s\n' "${FAILED[@]}"
    exit 1
fi
say "probe 10: all checks passed"
