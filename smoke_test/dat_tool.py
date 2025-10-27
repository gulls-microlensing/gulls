#!/usr/bin/env python3
"""Utility for inspecting and editing whitespace-delimited .dat catalog files."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple


class DatFileError(RuntimeError):
    """Raised when a .dat file cannot be parsed or written."""


def read_dat_file(path: Path) -> Tuple[str, List[str], List[List[str]]]:
    """Return header line, column names, and rows from a .dat file."""
    raw_lines = path.read_text().splitlines()

    # Drop leading blank lines to simplify parsing.
    while raw_lines and not raw_lines[0].strip():
        raw_lines.pop(0)

    if not raw_lines:
        raise DatFileError(f"{path} is empty")

    header_line = raw_lines[0]
    columns = header_line.split()
    if not columns:
        raise DatFileError(f"{path} has an empty header row")

    rows: List[List[str]] = []
    for line_number, line in enumerate(raw_lines[1:], start=2):
        if not line.strip():
            continue  # allow blank lines inside the file
        parts = line.split()
        if len(parts) != len(columns):
            raise DatFileError(
                f"{path}:{line_number} has {len(parts)} columns; "
                f"expected {len(columns)} based on header"
            )
        rows.append(parts)

    return header_line, columns, rows


def write_dat_file(path: Path, header: str, rows: Sequence[Sequence[str]]) -> None:
    """Persist a .dat file with the provided header and rows."""
    lines = [header, *((" ".join(row)) for row in rows)]
    data = "\n".join(lines) + "\n"
    path.write_text(data)


def resolve_column(columns: Sequence[str], index: int | None, name: str | None) -> Tuple[int, str]:
    """Return the column index and name from either an index or column label."""
    if index is not None and name is not None:
        raise DatFileError("Provide only one of --index or --column")
    if index is None and name is None:
        raise DatFileError("Either --index or --column must be specified")

    if index is not None:
        if index < 0 or index >= len(columns):
            raise DatFileError(f"Column index {index} is out of range (0-{len(columns) - 1})")
        col_name = columns[index]
        return index, col_name

    assert name is not None  # for type-checkers
    if name not in columns:
        suggestions = ", ".join(columns)
        raise DatFileError(f"Column '{name}' not found. Available columns: {suggestions}")
    resolved_index = columns.index(name)
    return resolved_index, name


def format_value(template: str, value: float, original: str) -> str:
    """Apply formatting template to a numeric value."""
    return template.format(value=value, original=original)


def perform_operation(existing: float, op: str, operand: float) -> float:
    """Compute the new value for a column."""
    if op in {"scale", "multiply", "mul"}:
        return existing * operand
    if op in {"add", "plus"}:
        return existing + operand
    if op == "set":
        return operand
    if op in {"pow", "power"}:
        return existing**operand
    raise DatFileError(f"Unsupported operation '{op}'. Choose from: multiply, add, set, power.")


def cmd_show(args: argparse.Namespace) -> int:
    header_line, columns, rows = read_dat_file(args.file)
    sample_row: List[str] | None = None

    if args.no_values:
        sample_row = None
    elif args.row is not None:
        if args.row < 0 or args.row >= len(rows):
            raise DatFileError(
                f"Row {args.row} is out of range. Available data rows: 0-{max(len(rows) - 1, 0)}"
            )
        sample_row = rows[args.row]
    elif rows:
        sample_row = rows[0]

    index_width = len(str(len(columns) - 1))
    name_width = max(len(name) for name in columns)

    print(f"{'Idx':>{index_width}}  {'Column':<{name_width}}", end="")
    if sample_row is not None:
        print("  SampleValue")
    else:
        print()

    for idx, name in enumerate(columns):
        print(f"{idx:>{index_width}}  {name:<{name_width}}", end="")
        if sample_row is not None:
            value = sample_row[idx]
            print(f"  {value}")
        else:
            print()

    if args.show_header:
        print("\nHeader:")
        print(header_line)

    return 0


def cmd_apply(args: argparse.Namespace) -> int:
    header_line, columns, rows = read_dat_file(args.file)
    col_index, col_name = resolve_column(columns, args.index, args.column)

    if not rows:
        raise DatFileError(f"{args.file} does not contain any data rows to modify")

    updated_rows: List[List[str]] = []
    changes: List[Tuple[int, str, str]] = []

    for row_number, current_row in enumerate(rows):
        current_value_str = current_row[col_index]
        current_value = float(current_value_str)

        new_value_number = perform_operation(current_value, args.operation, args.value)
        new_value_str = format_value(args.format, new_value_number, current_value_str)

        if args.precision_guard is not None:
            diff = abs(new_value_number - current_value)
            if diff < args.precision_guard:
                new_value_str = current_value_str

        new_row = list(current_row)
        new_row[col_index] = new_value_str
        updated_rows.append(new_row)

        if len(changes) < args.preview:
            changes.append((row_number, current_value_str, new_value_str))

    if args.inplace and not args.output:
        destination = args.file
    elif args.output is not None:
        destination = args.output
    else:
        raise DatFileError("Provide --output or use --inplace to write changes")

    write_dat_file(destination, header_line, updated_rows)

    target_desc = f"column '{col_name}' (index {col_index})"
    print(f"Applied {args.operation} to {target_desc} for {len(updated_rows)} rows.")
    print(f"Output written to {destination}")

    if changes:
        print("\nPreview of first modified rows:")
        for row_number, old, new in changes:
            print(f"  row {row_number}: {old} -> {new}")

    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Inspect and manipulate whitespace-delimited lens/source catalog files."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    show_parser = subparsers.add_parser(
        "show", help="List the columns and optional sample values for a .dat file."
    )
    show_parser.add_argument("file", type=Path, help="Path to the .dat file to inspect")
    show_parser.add_argument(
        "--row",
        type=int,
        default=None,
        help="Display values from the specified data row (0-indexed). Defaults to the first row.",
    )
    show_parser.add_argument(
        "--no-values",
        action="store_true",
        help="Only show column names without sample values.",
    )
    show_parser.add_argument(
        "--show-header",
        action="store_true",
        help="Print the raw header line after the column listing.",
    )
    show_parser.set_defaults(func=cmd_show)

    apply_parser = subparsers.add_parser(
        "apply", help="Modify a numeric column using simple arithmetic operations."
    )
    apply_parser.add_argument("file", type=Path, help="Path to the .dat file to modify")
    column_group = apply_parser.add_mutually_exclusive_group(required=True)
    column_group.add_argument("--column", type=str, help="Column name to operate on")
    column_group.add_argument("--index", type=int, help="Column index to operate on")
    apply_parser.add_argument(
        "--operation",
        type=str,
        choices=["add", "plus", "scale", "multiply", "mul", "set", "pow", "power"],
        required=True,
        help="Operation to apply to the column values.",
    )
    apply_parser.add_argument(
        "--value",
        type=float,
        required=True,
        help="Operand for the operation (ignored for operations that do not need one).",
    )
    apply_parser.add_argument(
        "--format",
        type=str,
        default="{value:.7e}",
        help="Python format string for the new values. You can reference '{value}' and '{original}'.",
    )
    apply_parser.add_argument(
        "--output",
        type=Path,
        help="Destination path for the modified file. Required unless --inplace is used.",
    )
    apply_parser.add_argument(
        "--inplace",
        action="store_true",
        help="Write changes directly back to the input file.",
    )
    apply_parser.add_argument(
        "--preview",
        type=int,
        default=5,
        help="Number of modified rows to preview after applying changes.",
    )
    apply_parser.add_argument(
        "--precision-guard",
        type=float,
        default=None,
        help="Skip rewriting values when the absolute change is below this threshold.",
    )
    apply_parser.set_defaults(func=cmd_apply)

    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
