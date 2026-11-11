#!/bin/sh
# End-to-end checks for the delaunayd binary. Usage: sh tests/cli_test.sh ./delaunayd
set -e
BIN=${1:-./delaunayd}
DIR=$(cd "$(dirname "$0")/.." && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
fail() { echo "FAIL: $*"; exit 1; }
pass=0
ok() { pass=$((pass + 1)); }

# 6x6 grid: 2 * 5 * 5 triangles
"$BIN" triangulate "$DIR/examples/grid.txt" > "$TMP/grid.txt"
head -1 "$TMP/grid.txt" | grep -q "36 points (36 distinct), 50 triangles, 85 edges" || fail "grid counts"; ok

# stdin input and -o output
"$BIN" triangulate -o "$TMP/out.txt" < "$DIR/examples/grid.txt"
cmp -s "$TMP/grid.txt" "$TMP/out.txt" || fail "stdin/-o mismatch"; ok

# constrained star: all 10 outline segments present
"$BIN" triangulate "$DIR/examples/star.txt" | head -1 | grep -q "10 constrained" || fail "star constraints"; ok

# collinear input: no triangles, strip-shaped Voronoi cells (duplicate ignored)
"$BIN" triangulate "$DIR/examples/collinear.txt" | head -1 | grep -q "6 points (5 distinct), 0 triangles, 4 edges" || fail "collinear"; ok
"$BIN" voronoi "$DIR/examples/collinear.txt" | grep -c '^line ' | grep -qx 4 || fail "collinear voronoi"; ok

# random points round-trip through the parser; Euler relation for the hull
"$BIN" random 500 --seed 9 > "$TMP/r.txt"
[ "$(grep -vc '^#' "$TMP/r.txt")" -eq 500 ] || fail "random count"; ok
H=$("$BIN" hull "$TMP/r.txt" | head -1 | sed 's/.*hull: \([0-9]*\) vertices.*/\1/')
T=$("$BIN" triangulate "$TMP/r.txt" | head -1 | sed 's/.*distinct), \([0-9]*\) triangles.*/\1/')
[ "$T" -eq $((2 * 500 - 2 - H)) ] || fail "euler: T=$T H=$H"; ok

# voronoi: one cell per site
"$BIN" voronoi "$TMP/r.txt" | grep -c '^site ' | grep -qx 500 || fail "voronoi cells"; ok

# locate
printf '0 0\n4 0\n0 4\n' > "$TMP/tri.txt"
"$BIN" locate 1 1 "$TMP/tri.txt" | grep -Eq "^triangle 0 (0 1 2|1 2 0|2 0 1)$" || fail "locate inside"; ok
[ "$("$BIN" locate 2 0 "$TMP/tri.txt")" = "edge 0 1" ] || fail "locate edge"; ok
[ "$("$BIN" locate 9 9 "$TMP/tri.txt")" = "outside" ] || fail "locate outside"; ok

# SVG output
for cmd in triangulate hull voronoi; do
  "$BIN" $cmd --svg "$DIR/examples/star.txt" > "$TMP/$cmd.svg" 2>/dev/null || [ $cmd = voronoi ]
done
grep -q '^<svg' "$TMP/triangulate.svg" || fail "svg triangulate"; ok
grep -q '^<svg' "$TMP/hull.svg" || fail "svg hull"; ok
"$BIN" voronoi --svg --delaunay "$TMP/r.txt" | grep -q 'stroke-dasharray' || fail "svg voronoi"; ok

# errors
if "$BIN" frobnicate 2>/dev/null; then fail "unknown command accepted"; fi; ok
printf '0 0\n1 1\n1 0\n0 1\nseg 0 1\nseg 2 3\n' | { if "$BIN" triangulate 2>"$TMP/err"; then fail "crossing accepted"; fi; }
grep -q crosses "$TMP/err" || fail "crossing message"; ok
printf '0 0\nfoo bar\n' | { if "$BIN" hull 2>"$TMP/err"; then fail "bad input accepted"; fi; }
grep -q 'line 2' "$TMP/err" || fail "parse error line"; ok

echo "cli tests passed ($pass checks)"
