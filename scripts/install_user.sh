#!/usr/bin/env bash
# Installation utilisateur hors-ligne, sans droits administrateur ni pip.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="$HOME/.local/bin"
PYZ="$ROOT/alignview.pyz"
[[ -r "$PYZ" ]] || { echo "Fichier manquant : $PYZ" >&2; exit 1; }
command -v python3 >/dev/null || { echo 'Python 3.10+ nécessaire' >&2; exit 1; }
python3 - <<'PY'
import sys
if sys.version_info < (3,10):raise SystemExit('AlignView nécessite Python >= 3.10')
PY
mkdir -p "$TARGET"
if [[ -e "$TARGET/alignview" ]] && [[ ${1:-} != --force ]]; then
  echo "Lancement déjà installé : $TARGET/alignview ; ajouter --force pour remplacer." >&2; exit 2
fi
cp "$PYZ" "$TARGET/alignview"
chmod 755 "$TARGET/alignview"
echo "Installé dans $TARGET/alignview"
echo "Ajouter $TARGET au PATH si nécessaire : export PATH=\"\$HOME/.local/bin:\$PATH\""
echo 'Lancer : alignview --version'
