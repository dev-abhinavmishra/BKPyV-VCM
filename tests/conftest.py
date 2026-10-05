"""Test import-path setup for the repository and local utility scripts."""

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "scripts"
SRC_DIR = REPO_ROOT / "src"

for path in (REPO_ROOT, SCRIPTS_DIR, SRC_DIR):
    path_string = str(path)
    if path_string not in sys.path:
        sys.path.insert(0, path_string)
