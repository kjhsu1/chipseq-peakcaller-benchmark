"""Pytest path setup for local helper-module imports."""

"""Imports"""

import sys
from pathlib import Path


"""Module Setup"""

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
