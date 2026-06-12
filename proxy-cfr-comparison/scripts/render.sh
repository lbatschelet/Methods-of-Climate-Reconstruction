#!/usr/bin/env bash
# Render analysis.qmd with the project virtualenv (overrides broken global QUARTO_PYTHON).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# Override a broken global QUARTO_PYTHON (e.g. ~/.quarto-env/bin/python)
unset QUARTO_PYTHON
export QUARTO_PYTHON="$ROOT/.venv/bin/python"
cd "$ROOT"
bash scripts/ensure_quarto_python.sh >/dev/null
exec quarto render "$@"
