#!/usr/bin/env bash
# Create .venv, install Python deps, register Jupyter kernel for Quarto.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

PY="${PYTHON:-}"
if [[ -z "$PY" ]]; then
  for candidate in \
    /opt/homebrew/Cellar/python@3.14/*/Frameworks/Python.framework/Versions/3.14/bin/python3.14 \
    /opt/homebrew/bin/python3 \
    /usr/bin/python3; do
    if [[ -x "$candidate" ]] && "$candidate" -c "import sys; print(sys.version)" >/dev/null 2>&1; then
      PY="$candidate"
      break
    fi
  done
fi
if [[ -z "$PY" ]] || ! command -v "$PY" >/dev/null 2>&1; then
  echo "Error: no working Python 3 found. Install Python 3.11+ or set PYTHON=..." >&2
  exit 1
fi
if ! "$PY" -c "import sys; print(sys.version)" >/dev/null 2>&1; then
  echo "Error: $PY is broken (possibly overwritten). Try: brew reinstall python@3.14" >&2
  exit 1
fi

echo "==> Using Python: $PY"

echo "==> Creating virtualenv at $ROOT/.venv"
"$PY" -m venv .venv
.venv/Scripts/pip install --upgrade pip
.venv/Scripts/pip install -r python/requirements.txt

SITE_PACKAGES="$(.venv/bin/python -c "import site; print(site.getsitepackages()[0])")"
echo "$ROOT/python" > "$SITE_PACKAGES/proxy_cfr_project.pth"
echo "==> Added $ROOT/python to venv path (.pth)"

echo "==> Registering Jupyter kernel 'dod2k-cfr'"
.venv/bin/python -m ipykernel install --user --name=dod2k-cfr --display-name="Python (dod2k-cfr)"

bash "$ROOT/scripts/ensure_quarto_python.sh"

echo ""
echo "Done. Plain 'quarto render' should work from proxy-cfr-comparison/."
echo "Optional: clone dod2k utilities — bash scripts/setup_dod2k_repo.sh"
echo ""
echo "R package for Parquet bridge:"
echo "  Rscript -e \"install.packages('arrow', repos='https://cloud.r-project.org')\""
