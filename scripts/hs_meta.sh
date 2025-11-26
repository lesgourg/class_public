set -euo pipefail
gitshort="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
now="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
ccver="$(${CC:-cc} --version 2>/dev/null | head -n1 || echo 'cc: unknown')"
pycmd="${PYTHON:-python3}"
if command -v "$pycmd" >/dev/null 2>&1; then
  pyver="$("$pycmd" --version 2>&1)"
else
  pyver="python: unknown"
fi
printf $'rev: %s\ndate: %s\ncc: %s\npython: %s\n' "$gitshort" "$now" "$ccver" "$pyver" > hs_meta.txt
