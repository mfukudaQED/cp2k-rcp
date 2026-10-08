#!/usr/bin/env python3
"""Validate native CP2K RCP k-point output and symmetry consistency."""

import argparse
import math
from pathlib import Path


def metric(path: Path, label: str) -> float:
    """Return the numeric value of an RCP diagnostic in a CP2K output."""
    for line in path.read_text().splitlines():
        if line.startswith(" RCP|") and label in line:
            return float(line.split()[-1])
    raise ValueError(f"Missing metric {label} in {path}")


def read_cube(path: Path) -> list[float]:
    """Read Gaussian cube data and reject malformed or non-finite values."""
    lines = path.read_text().splitlines()
    natoms = abs(int(lines[2].split()[0]))
    numbers = []
    for line in lines[6 + natoms :]:
        for item in line.split():
            number = float(item.replace("D", "E"))
            if not math.isfinite(number):
                raise ValueError(f"Nonfinite cube value in {path}")
            numbers.append(number)
    expected = math.prod(abs(int(lines[i].split()[0])) for i in range(3, 6))
    if len(numbers) != expected:
        raise ValueError(f"Wrong cube point count in {path}: {len(numbers)} vs {expected}")
    return numbers


def validate(root: Path) -> None:
    """Check RCP charge and energy diagnostics across four k-point cases."""
    cases = (
        "H2-rcp-kpoints",
        "H2-rcp-kpoints-sym",
        "H2-rcp-kpoints-1x1x1",
        "H2-rcp-kpoints-uks",
        "H2-rcp-kpoints-parallel",
        "H2-rcp-gamma-periodic",
    )
    for name in cases:
        out = root / f"{name}.out"
        content = out.read_text()
        if "PROGRAM ENDED AT" not in content:
            raise AssertionError(f"CP2K did not finish: {out}")
        a = metric(out, "Electron count from sum(weights)")
        b = metric(out, "Electron count from Tr[P_window*S]")
        c = metric(out, "Electron count from grid integration")
        kinetic_error = metric(out, "|Integral(T_e)-E_kin|")
        gradient_error = metric(out, "|Laplacian-gradient integrals|")
        assert abs(a - 2.0) < 1e-6, (name, a)
        assert abs(a - b) < 1e-6, (name, a, b)
        assert abs(b - c) < 1e-6, (name, b, c)
        assert kinetic_error < 1e-4, (name, kinetic_error)
        assert gradient_error < 1e-4, (name, gradient_error)
        expected_spin = 2 if name.endswith("uks") else 1
        assert metric(out, "Number of collinear spin channels") == expected_spin
        print(f"PASS {name}: Q={c:.12f}, kinetic error={kinetic_error:.3e}")

    for suffix in (
        "RCP1_0.cube",
        "RCP-KINETIC-ENERGY-DENSITY1_0.cube",
        "RCP-REGIONAL-ENERGY-DENSITY1_0.cube",
    ):
        full = read_cube(root / f"H2-rcp-kpoints-{suffix}")
        reduced = read_cube(root / f"H2-rcp-kpoints-sym-{suffix}")
        if len(full) != len(reduced):
            raise AssertionError(f"Different cube grids: {suffix}")
        max_diff = max(abs(a - b) for a, b in zip(full, reduced))
        assert max_diff < 1e-7, (suffix, max_diff)
        print(f"PASS symmetry cube {suffix}: max absolute difference={max_diff:.3e}")

    # Check the selected-density cube from the inter-k-point MPI-group test.
    density_values = read_cube(root / "H2-rcp-kpoints-parallel-RCP-DENSITY-WINDOW1_0.cube")
    if min(density_values) < -1e-8:
        raise AssertionError("Window electron density contains negative values")
    print(f"PASS k-point MPI-group density cube: {len(density_values)} finite grid values")

    # The Gamma-point and one-k-point branches should agree on the same periodic cell.
    gamma = root / "H2-rcp-gamma-periodic.out"
    one_k = root / "H2-rcp-kpoints-1x1x1.out"
    e_gamma = metric(gamma, "Regional energy density integral")
    e_one_k = metric(one_k, "Regional energy density integral")
    assert abs(e_gamma - e_one_k) < 1e-6, (e_gamma, e_one_k)
    print(f"PASS Gamma/1x1x1 regional energy consistency: delta={abs(e_gamma-e_one_k):.3e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Directory with CP2K outputs and cubes")
    args = parser.parse_args()
    validate(args.directory)
