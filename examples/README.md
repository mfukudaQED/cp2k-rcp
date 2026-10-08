# CP2K RCP example input decks

Every `*.inp` file here is intended for **CP2K 2026.2 with the RCP patch installed**. These directories contain text inputs only. No computed wavefunctions, restart files, CUBE output, or simulation logs are shipped.

## Directory guide

| Directory | Input(s) | Purpose |
| --- | --- | --- |
| `h2/` | `H2-rcp.inp` | Isolated H₂ molecule, Gamma-point, closed-shell RKS |
| `h2/` | `H2-rcp-gamma-periodic.inp` | Periodic H₂ reference with the Gamma-only execution path |
| `h2/` | `H2-rcp-kpoints-1x1x1.inp` | Same periodic cell with an explicit 1×1×1 k-point path |
| `h2/` | `H2-rcp-kpoints.inp` | Monkhorst–Pack 2×1×1, full grid, complex orbitals |
| `h2/` | `H2-rcp-kpoints-sym.inp` | Symmetry-reduced 2×1×1 grid |
| `h2/` | `H2-rcp-kpoints-parallel.inp` | k-point MPI-group path |
| `h2/` | `H2-rcp-kpoints-uks.inp` | Collinear UKS with multiplicity 1 (code-path regression, **not** an open-shell example) |
| `benzene/` | `benzene.inp` | Isolated benzene molecule, RKS, −3 to 0 eV RCP window |
| `c2h5/` | `c2h5.inp` | Isolated C₂H₅ radical, UKS doublet, −2 to 0 eV RCP window |

The H₂ reference expectations are collected in `h2/TEST_FILES.toml`. For the real open-shell case, use the C₂H₅ input, not the artificial UKS-singlet regression.

## Running an example

Build a **patched** CP2K 2026.2 binary and ensure its standard basis and potential data files are installed. In a terminal in the directory containing the input:

~~~bash
cd h2
cp2k.psmp -i H2-rcp.inp -o H2-rcp.out
grep 'RCP|' H2-rcp.out
ls -lh *RCP*.cube
~~~

The H₂ examples default to `PRINT_DENSITY_WINDOW F`, so usually three CUBE types are written. Set that keyword to `T` to include the energy-window density CUBE as well.

For an MPI Slurm job on sham:

~~~bash
cd h2
sbatch --export=ALL,CP2K_ENV=/path/to/patched-cp2k/install/cp2k_env \
  ../../tools/run_cp2k_sham.slurm H2-rcp.inp H2-rcp.out
~~~

Adjust job resources in the Slurm template to suit the calculation.

## K-point regression data layout

`tools/check_kpoint_consistency.py` validates *computed* results. It expects the named k-point cases to have been executed with output files named after their input stems, e.g.:

~~~text
regression-results/
  H2-rcp-kpoints.out
  H2-rcp-kpoints-RCP1_0.cube
  H2-rcp-kpoints-RCP-KINETIC-ENERGY-DENSITY1_0.cube
  H2-rcp-kpoints-RCP-REGIONAL-ENERGY-DENSITY1_0.cube
  ...
~~~

The validator expects all six periodic reference cases (including symmetry, UKS, and MPI-group cases). Place their outputs in **one directory** named `regression-results/`, then run:

~~~bash
python3 tools/check_kpoint_consistency.py regression-results/
~~~

This directory is deliberately **not provided** and should stay outside Git.

## Limitations

- These are small demonstrations and regression inputs, not converged predictive materials calculations.
- RCP CUBE values have physical units and require a common-grid comparison and a mask in low electron-density regions.
- Upstream CP2K without the RCP patch cannot parse `&DFT / &PRINT / &RCP`.
- A larger 210-carbon restart-dependent chain and a reduced diamond(001) model exist in the internal development workspace, but are **not included** in this public example bundle. Their structures, provenance, and associated reference outputs require a separate release decision.
