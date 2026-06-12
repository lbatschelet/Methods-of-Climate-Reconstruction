#!/usr/bin/env bash
# Create .venv, install Python deps, register Jupyter kernel for Quarto.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

PY="${PYTHON:-python3}"
if ! command -v "$PY" >/dev/null 2>&1; then
  echo "Error: $PY not found. Install Python 3.11+ or set PYTHON=..." >&2
  exit 1
fi

echo "==> Creating virtualenv at $ROOT/.venv"
"$PY" -m venv .venv
.venv/bin/pip install --upgrade pip
.venv/bin/pip install -r python/requirements.txt

SITE_PACKAGES="$(.venv/bin/python -c "import site; print(site.getsitepackages()[0])")"
echo "$ROOT/python" > "$SITE_PACKAGES/proxy_cfr_project.pth"
echo "==> Added $ROOT/python to venv path (.pth)"

echo "==> Registering Jupyter kernel 'dod2k-cfr'"
.venv/bin/python -m ipykernel install --user --name=dod2k-cfr --display-name="Python (dod2k-cfr)"

echo "QUARTO_PYTHON=$ROOT/.venv/bin/python" > "$ROOT/_environment"

echo ""
echo "Done. _environment points QUARTO_PYTHON to .venv/bin/python."
echo "Optional: clone dod2k utilities — bash scripts/setup_dod2k_repo.sh"
echo ""
echo "R package for Parquet bridge:"
echo "  Rscript -e \"install.packages('arrow', repos='https://cloud.r-project.org')\""
