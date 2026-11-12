#!/usr/bin/env bash
# Offline user installation without administrator access or pip.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="$HOME/.local/bin"
PYZ="$ROOT/alignview.pyz"
[[ -r "$PYZ" ]] || { echo "Missing file : $PYZ" >&2; exit 1; }
command -v python3 >/dev/null || { echo 'Python 3.10+ is required' >&2; exit 1; }
python3 - <<'PY'
import sys
if sys.version_info < (3,10):raise SystemExit('AlignView requires Python >= 3.10')
PY
mkdir -p "$TARGET"
if [[ -e "$TARGET/alignview" ]] && [[ ${1:-} != --force ]]; then
  echo "Launcher already installed : $TARGET/alignview ; add --force to replace it." >&2; exit 2
fi
cp "$PYZ" "$TARGET/alignview"
chmod 755 "$TARGET/alignview"
echo "Installed to $TARGET/alignview"
echo "Add $TARGET to PATH if needed : export PATH=\"\$HOME/.local/bin:\$PATH\""
echo 'Run: alignview --version'
