#!/usr/bin/env bash
# Clone the official DoD2k repo for dod2k_utilities (optional but recommended).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TARGET="${ROOT}/external/dod2k"
if [[ -d "${TARGET}/dod2k_utilities" ]]; then
  echo "Already present: ${TARGET}"
  exit 0
fi
mkdir -p "${ROOT}/external"
git clone --depth 1 https://github.com/lluecke/dod2k.git "${TARGET}"
echo "Cloned to ${TARGET}"
echo "Set: export DOD2K_REPO=${TARGET}"
