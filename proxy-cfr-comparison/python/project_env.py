"""Resolve project paths for Quarto {python} chunks."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def find_project_root() -> Path:
    env = os.environ.get("PROXY_CFR_PROJECT_ROOT")
    if env:
        p = Path(env).resolve()
        if (p / "config" / "experiments.R").exists():
            return p

    candidates = [Path.cwd()]
    if (candidates[0] / "analysis.qmd").exists():
        return candidates[0].resolve()
    if (candidates[0] / "proxy-cfr-comparison" / "config" / "experiments.R").exists():
        return (candidates[0] / "proxy-cfr-comparison").resolve()

    for parent in [Path.cwd(), *Path.cwd().parents]:
        for root in (parent, parent / "proxy-cfr-comparison"):
            if (root / "config" / "experiments.R").exists():
                return root.resolve()

    raise FileNotFoundError(
        "Cannot find proxy-cfr-comparison/. Run Quarto from the project directory."
    )


def bootstrap_python_path(root: Path | None = None) -> Path:
    root = root or find_project_root()
    py_dir = root / "python"
    if str(py_dir) not in sys.path:
        sys.path.insert(0, str(py_dir))

    repo = os.environ.get(
        "DOD2K_REPO",
        str(root.parent / "external" / "dod2k"),
    )
    if Path(repo).joinpath("dod2k_utilities").is_dir() and repo not in sys.path:
        sys.path.insert(0, repo)

    return root


def data_paths(root: Path | None = None) -> dict[str, Path]:
    root = root or find_project_root()
    repo = root.parent
    return {
        "project": root,
        "repo": repo,
        "dod2k_csv": repo / "course" / "data" / "dod2k_v2.0",
        "cache": root / "output" / "cache",
    }
