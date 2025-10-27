#!/usr/bin/env python3
"""Assemble the GitHub release body from optional release notes and plot notes."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

DEFAULT_BODY = textwrap.dedent(
    """\
    ## ⚠️ Important Notice

    **Binary releases include GSL fallback implementations and are NOT suitable for production science.**

    - **For production use**: Download source code and compile with licensed Numerical Recipes
    - **For testing/development**: Binary releases are fine for testing workflows
    - **See documentation**: Details on replacing GSL stubs with Numerical Recipes

    ## What's Included

    - **Source code**: Complete Gulls source with CMake build system
    - **Binaries**: Linux executables (GSL fallbacks - testing only)
    - **Documentation**: Built HTML documentation
    - **Smoke test plots**: Visual proof that the release works (see below)
    """
)


def _read_section(path: Path) -> str:
    if path.is_file():
        content = path.read_text(encoding="utf-8").strip()
        if content:
            return content
    return ""


def build_release_body(release_notes: Path, plot_notes: Path) -> str:
    sections = []

    release_section = _read_section(release_notes)
    if release_section:
        sections.append(release_section)
    else:
        sections.append(DEFAULT_BODY.strip())

    plot_section = _read_section(plot_notes)
    if plot_section:
        sections.append(plot_section)

    return "\n\n".join(sections).strip() + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--release-notes",
        default="RELEASE_NOTES.md",
        help="Path to the base release notes (default: %(default)s)",
    )
    parser.add_argument(
        "--plots",
        default="PLOT_NOTES.md",
        help="Path to the generated plot notes (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        default="RELEASE_BODY.md",
        help="Path to write the final release body (default: %(default)s)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    release_notes_path = Path(args.release_notes)
    plot_notes_path = Path(args.plots)
    output_path = Path(args.output)

    body = build_release_body(release_notes_path, plot_notes_path)
    output_path.write_text(body, encoding="utf-8")


if __name__ == "__main__":
    main()
