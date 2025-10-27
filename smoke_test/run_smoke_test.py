#!/usr/bin/env python3
"""Minimal entry-point wrapper for the gulls smoke test."""
from __future__ import annotations

import sys
from pathlib import Path

if __package__ is None or __package__ == "":
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from smoke_test.runner import main  # pylint: disable=wrong-import-position


if __name__ == "__main__":
    sys.exit(main())
