#!/usr/bin/env bash
# Ensure Quarto uses this project's .venv (fixes broken global QUARTO_PYTHON).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VENV_PY="$ROOT/.venv/bin/python"

if [[ ! -x "$VENV_PY" ]]; then
  echo "==> No .venv found — running setup_python.sh"
  bash "$ROOT/scripts/setup_python.sh"
fi

# Global Quarto installs often set QUARTO_PYTHON=~/.quarto-env/bin/python.
# Write a small wrapper (never touch .venv/bin/python or Homebrew Python).
QUARTO_ENV="$HOME/.quarto-env/bin/python"
mkdir -p "$(dirname "$QUARTO_ENV")"
rm -f "$QUARTO_ENV"
cat > "$QUARTO_ENV" <<EOF
#!/usr/bin/env bash
exec "$VENV_PY" "\$@"
EOF
chmod +x "$QUARTO_ENV"
echo "==> Updated $QUARTO_ENV wrapper -> $VENV_PY"

"$VENV_PY" -m ipykernel install --user --name=dod2k-cfr --display-name="Python (dod2k-cfr)" >/dev/null

echo "QUARTO_PYTHON=$VENV_PY" > "$ROOT/_environment"
