#!/usr/bin/env bash
set -euo pipefail
ALIGNVIEW_SCRIPT="${BASH_SOURCE[0]}"
while [[ -L "$ALIGNVIEW_SCRIPT" ]]; do
  ALIGNVIEW_LINK_DIR="$(cd -- "$(dirname -- "$ALIGNVIEW_SCRIPT")" && pwd)"
  ALIGNVIEW_SCRIPT="$(readlink -- "$ALIGNVIEW_SCRIPT")"
  [[ "$ALIGNVIEW_SCRIPT" == /* ]] || ALIGNVIEW_SCRIPT="$ALIGNVIEW_LINK_DIR/$ALIGNVIEW_SCRIPT"
done
ALIGNVIEW_ROOT="$(cd -- "$(dirname -- "$ALIGNVIEW_SCRIPT")" && pwd)"
if [[ ${1:-} == --classic ]]; then
  shift
  [[ -f "$ALIGNVIEW_ROOT/archives/alignview_legacy.sh" ]] || {
    echo 'Interface historique absente de cette distribution. Lancer AlignView sans --classic.' >&2; exit 2;
  }
  exec bash "$ALIGNVIEW_ROOT/archives/alignview_legacy.sh" "$@"
fi
if ! command -v python3 >/dev/null; then echo 'Python 3.10+ requis' >&2; exit 127; fi
# Detect an optional local Conda prefix. An explicit
# environment variable or --tools-prefix still takes priority.
if [[ -z ${ALIGNVIEW_TOOLS_PREFIX:-} && -d "$ALIGNVIEW_ROOT/tools/conda-meta" ]]; then
  export ALIGNVIEW_TOOLS_PREFIX="$ALIGNVIEW_ROOT/tools"
fi
# Keep the caller's directory: relative inputs and outputs belong to that place.
exec python3 "$ALIGNVIEW_ROOT/alignview.py" "$@"
