#!/usr/bin/env python3
"""Normalize Gaussian cube volumetric data so every value is whitespace separated.

Legacy CP2K cube output used fixed-width E13.5E3 fields with no explicit
separator. Consecutive negative values therefore touched each other and some
viewers could not parse the file. This utility preserves every numeric token
exactly and only rewrites the volumetric section using one explicit leading
space plus a right-justified 13-character value, with at most six values per
line.

Only scalar cube files (non-negative atom count) are modified. Files are
replaced atomically after point-count and token-hash verification.
"""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import re
import stat
import sys
import tempfile
from typing import Iterable

NUMBER_RE = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?")
VALUES_PER_LINE = 6
FIELD_WIDTH = 13


def iter_cube_files(paths: Iterable[Path]) -> list[Path]:
    """Return unique real cube files below paths, skipping symlinks."""
    found: set[Path] = set()
    for path in paths:
        path = path.expanduser()
        if path.is_symlink():
            continue
        if path.is_file():
            if path.suffix.lower() == ".cube":
                found.add(path.resolve())
            continue
        if path.is_dir():
            for candidate in path.rglob("*.cube"):
                if candidate.is_symlink() or not candidate.is_file():
                    continue
                found.add(candidate.resolve())
    return sorted(found)


def token_digest_update(digest, token: str) -> None:
    """Update a whitespace-independent token digest."""
    digest.update(token.encode("ascii"))
    digest.update(b"\n")


def parse_header(src) -> tuple[list[str], int]:
    """Read a scalar-cube header and return it with expected grid point count."""
    header: list[str] = []

    for _ in range(2):
        line = src.readline()
        if not line:
            raise ValueError("truncated cube header")
        header.append(line)

    atoms_line = src.readline()
    if not atoms_line:
        raise ValueError("missing atom-count line")
    header.append(atoms_line)
    parts = atoms_line.split()
    if not parts:
        raise ValueError("invalid atom-count line")
    natoms = int(parts[0])
    if natoms < 0:
        raise ValueError("multi-dataset/orbital cubes with negative atom count are not supported")

    shape: list[int] = []
    for _ in range(3):
        line = src.readline()
        if not line:
            raise ValueError("truncated grid header")
        header.append(line)
        parts = line.split()
        if not parts:
            raise ValueError("invalid grid line")
        shape.append(abs(int(parts[0])))

    for _ in range(natoms):
        line = src.readline()
        if not line:
            raise ValueError("truncated atom list")
        header.append(line)

    return header, shape[0] * shape[1] * shape[2]


def scan_cube(path: Path) -> tuple[int, str, bool]:
    """Return point count, token digest and whether adjacent tokens have whitespace."""
    digest = hashlib.sha256()
    count = 0
    spacing_ok = True

    with path.open("r", encoding="ascii", errors="strict") as src:
        _, expected = parse_header(src)
        for line in src:
            matches = list(NUMBER_RE.finditer(line))
            for idx, match in enumerate(matches):
                if idx > 0:
                    between = line[matches[idx - 1].end() : match.start()]
                    if not between or not any(ch.isspace() for ch in between):
                        spacing_ok = False
                token = match.group(0)
                token_digest_update(digest, token)
                count += 1

    if count != expected:
        raise ValueError(f"volumetric point count {count} != expected {expected}")
    return count, digest.hexdigest(), spacing_ok


def normalize_cube(path: Path, *, dry_run: bool = False) -> tuple[int, int, int, bool]:
    """Rewrite one cube atomically and verify that all numeric tokens are unchanged."""
    old_size = path.stat().st_size
    point_count, original_hash, already_spacing_ok = scan_cube(path)

    if dry_run:
        return point_count, old_size, old_size, already_spacing_ok

    mode = stat.S_IMODE(path.stat().st_mode)
    tmp_path: Path | None = None

    try:
        with path.open("r", encoding="ascii", errors="strict") as src:
            header, expected = parse_header(src)
            fd, tmp_name = tempfile.mkstemp(
                prefix=f".{path.name}.spacing.",
                suffix=".tmp",
                dir=str(path.parent),
                text=True,
            )
            tmp_path = Path(tmp_name)

            with os.fdopen(fd, "w", encoding="ascii", newline="\n") as dst:
                dst.writelines(header)
                chunk: list[str] = []
                count = 0
                digest = hashlib.sha256()

                for line in src:
                    for match in NUMBER_RE.finditer(line):
                        token = match.group(0)
                        token_digest_update(digest, token)
                        count += 1
                        chunk.append(token)
                        if len(chunk) == VALUES_PER_LINE:
                            dst.write(
                                "".join(" " + item.rjust(FIELD_WIDTH) for item in chunk) + "\n"
                            )
                            chunk.clear()

                if chunk:
                    dst.write("".join(" " + item.rjust(FIELD_WIDTH) for item in chunk) + "\n")

                dst.flush()
                os.fsync(dst.fileno())

        if count != expected:
            raise ValueError(f"volumetric point count {count} != expected {expected}")
        if digest.hexdigest() != original_hash:
            raise ValueError("numeric token hash changed during conversion")

        os.chmod(tmp_path, mode)
        os.replace(tmp_path, path)
        tmp_path = None

        new_count, new_hash, spacing_ok = scan_cube(path)
        if new_count != point_count or new_hash != original_hash:
            raise ValueError("post-conversion numeric verification failed")
        if not spacing_ok:
            raise ValueError("post-conversion spacing verification failed")

        return point_count, old_size, path.stat().st_size, already_spacing_ok
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Add explicit whitespace separators to CP2K Gaussian cube data."
    )
    parser.add_argument("paths", nargs="+", type=Path, help="Cube file(s) or directories.")
    parser.add_argument(
        "--check",
        action="store_true",
        help="Only verify point counts and report spacing; do not modify files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show files that would be processed without modifying them.",
    )
    args = parser.parse_args()

    files = iter_cube_files(args.paths)
    if not files:
        print("No cube files found.", file=sys.stderr)
        return 1

    failures = 0
    for path in files:
        try:
            if args.check:
                points, _, spacing_ok = scan_cube(path)
                status = "OK" if spacing_ok else "NEEDS-SPACING"
                print(f"{status:13s} {points:10d}  {path}")
                continue

            points, old_size, new_size, was_ok = normalize_cube(path, dry_run=args.dry_run)
            status = "DRY-RUN" if args.dry_run else ("NORMALIZED" if not was_ok else "REWRITTEN")
            print(
                f"{status:10s} {points:10d} points  "
                f"{old_size / 1024**2:9.2f} -> {new_size / 1024**2:9.2f} MiB  {path}"
            )
        except Exception as exc:
            failures += 1
            print(f"ERROR {path}: {exc}", file=sys.stderr)

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
