#!/usr/bin/env bash
# Render analysis.qmd with the project virtualenv (overrides broken global QUARTO_PYTHON).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export QUARTO_PYTHON="$ROOT/.venv/bin/python"
cd "$ROOT"
exec quarto render analysis.qmd "$@"
