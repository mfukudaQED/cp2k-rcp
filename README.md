# Regional Chemical Potential (RCP) for CP2K

A research extension of **CP2K 2026.2** for calculating the regional chemical potential (RCP) and related real-space electron and kinetic-energy quantities from converged Kohn–Sham states.

This repository distributes a **patch against upstream CP2K**, user documentation, representative CP2K inputs, and small validation/utilities scripts. It deliberately does **not** mirror the CP2K source tree or host simulation output (such as CUBE volumetric files, wavefunction/restart files, or Slurm logs).

> **Status:** Research software, tested for GPW and collinear-spin electronic structure. It is not an official CP2K release or an upstream CP2K feature.

## Features and scope

| Feature | Status |
| --- | --- |
| Gaussian and plane waves (GPW) with GTH pseudopotentials | Supported; PBE extensively used in validation |
| Gamma-point calculations (molecules or periodic systems) | Supported |
| General k-point sampling with complex Bloch orbitals | Supported; tested for small meshes, symmetry reduction, and k-point MPI groups |
| Restricted Kohn–Sham (RKS) | Supported |
| Unrestricted, **collinear** Kohn–Sham (UKS) | Supported; spin channels are summed in the output fields |
| Energy-window selection relative to a reference chemical potential | Supported |
| Fermi–Dirac smearing in the parent SCF calculation | Supported; material-specific convergence studies still required |
| All-electron GAPW; ZORA/DKH scalar relativity | **Not implemented for RCP**; explicit runtime guard |
| Noncollinear spin, spin spirals, spin–orbit coupling | **Not implemented for RCP** |
| Hybrid/HF and general metallic k-mesh convergence | Not systematically validated |

The patch also contains an **isolated experimental Pauli-spinor local kernel** and unit tests. That module is **not integrated with CP2K's noncollinear/SOC self-consistent calculation or CUBE output** and should not be interpreted as support for those features.

### What is calculated?

For an energy window specified relative to the reference chemical potential, the extension outputs the energy-window electron density and regional energy density and their ratio:

```math
\mu_R^\tau(\mathbf r)
  = \frac{\varepsilon_{\tau,\mathrm{EW}}(\mathbf r)}
         {n_{\mathrm{EW}}(\mathbf r)}.
```

A separate output provides the Laplacian-form kinetic-energy density from **all occupied states**, useful, for example, for coloring an electronic isosurface with RCP values.

The density threshold masks poorly defined ratios in near-vacuum / orbital-node regions. The energy-window broadening is distinct from the parent SCF electronic smearing.

## Repository contents

- [patches/cp2k-2026.2-rcp.patch](patches/cp2k-2026.2-rcp.patch) — **complete patch**, recommended for installation.
- [patches/cp2k-2026.2-rcp-series.mbox](patches/cp2k-2026.2-rcp-series.mbox) — the same final changes as a **single squashed Git email patch**, installable with `git am` (it does not preserve the separate development commits). **Choose one installation method; do not apply both.**
- [docs/rcp_quickstart.md](docs/rcp_quickstart.md) — **start here if you are new to CP2K or theoretical/computational chemistry**: one-process H₂ example and explanations of output files.
- [docs/rcp_user_manual.md](docs/rcp_user_manual.md) — platform-independent user manual: executable setup, inputs, outputs, diagnostics, and validation.
- [docs/rcp_tutorial.md](docs/rcp_tutorial.md) — detailed RCP theory, numerical formulas, and historical validation notes.
- [docs/sham_notes.md](docs/sham_notes.md) — optional site-specific HPC notes, **not needed on other computers**.
- [examples/README.md](examples/README.md) — curated, fully specified CP2K inputs for H2, benzene, and the C2H5 radical.
- [tools/normalize_cube_spacing.py](tools/normalize_cube_spacing.py) — normalize spacing in older CP2K Gaussian CUBE files without repeating SCF calculations.
- [tools/check_kpoint_consistency.py](tools/check_kpoint_consistency.py) — check diagnostics and CUBE consistency **after** generating the expected H2 regression outputs.
- [tools/run_cp2k_sham.slurm](tools/run_cp2k_sham.slurm) — optional sham Slurm template without personal machine paths.
- [LICENSE](LICENSE) — GNU GPL version 2 license text; the CP2K-derived source code uses **GPL-2.0-or-later**.

No CUBE files or precomputed SCF results are distributed.

The patch also includes an **English in-tree copy of the RCP documentation** under CP2K's `docs/methods/`. The beginner guide and standalone manuals in this repository's `docs/`, along with the supplied `examples/`, are the preferred public entry points; the historical validation examples described in the detailed tutorial are not all distributed.

## Install the patch

The patch was generated from the `rcp-development` branch relative to:

- Upstream project: [cp2k/cp2k](https://github.com/cp2k/cp2k)
- CP2K release: **2026.2**
- Exact upstream base commit: **`67b5da876dd6a76b8b021d5a04d1c81ba79a4c50`**

Start with a clean CP2K checkout and a separate clone or download of this distribution:

~~~bash
git clone https://github.com/cp2k/cp2k.git cp2k
cd cp2k

git checkout -b rcp-patched 67b5da876dd6a76b8b021d5a04d1c81ba79a4c50
git apply --check ../cp2k-rcp/patches/cp2k-2026.2-rcp.patch
git apply ../cp2k-rcp/patches/cp2k-2026.2-rcp.patch
~~~

These commands assume you placed **this distribution** in `../cp2k-rcp`. Adjust the path if it is elsewhere. The `--check` step intentionally fails if you use an incompatible CP2K revision or a tree with conflicting changes.

Alternatively, from the same clean upstream commit:

~~~bash
git am ../cp2k-rcp/patches/cp2k-2026.2-rcp-series.mbox
~~~

**Do not use `git apply` and `git am` on the same checkout.**

Build the patched CP2K source using the [official CP2K build instructions](https://manual.cp2k.org/) appropriate to your environment. The patch is not a prebuilt executable.

## Run a minimal example

For a fully guided first calculation, see the [beginner's tutorial](docs/rcp_quickstart.md). After installing a patched CP2K executable:

~~~bash
cd /path/to/cp2k-rcp/examples/h2
cp2k.psmp -i H2-rcp.inp -o H2-rcp.out

grep 'SCF run converged' H2-rcp.out
grep 'RCP|' H2-rcp.out
ls *RCP*.cube
~~~

The example is an H2 **Gamma-point, RKS** single-point calculation. The H2 regression inputs intentionally disable the extra window-density CUBE by default to limit output volume; change `PRINT_DENSITY_WINDOW F` to `T` if you want that file.

For a periodic k-point test, replace `H2-rcp.inp` with `H2-rcp-kpoints.inp`. For an open-shell collinear UKS example, use `examples/c2h5/c2h5.inp`.

### Minimal RCP input block

Insert under `&FORCE_EVAL / &DFT / &PRINT`:

~~~text
&PRINT
  &RCP ON
    ENERGY_LOWER [eV] -3.0
    ENERGY_UPPER [eV] 0.0
    BROADENING_LOWER [eV] 0.001
    BROADENING_UPPER [eV] 0.001
    DENSITY_CUTOFF 1.0E-12

    PRINT_DENSITY_WINDOW T
    PRINT_KINETIC_ENERGY_DENSITY T
    PRINT_REGIONAL_ENERGY_DENSITY T
    PRINT_RCP T

    MPI_IO F
    STRIDE 1 1 1
  &END RCP
&END PRINT
~~~

This fragment is **not a standalone CP2K input**; use a complete file from `examples/`. The window bounds and broadening are in eV relative to the reference chemical potential. `STRIDE` affects CUBE output resolution. `MPI_IO F` selects a serial CUBE writer, **not** a serial SCF calculation.

## Output files

For `PROJECT mycalc`, representative CUBE filenames are:

| Filename pattern | Meaning | Data unit |
| --- | --- | --- |
| `mycalc-RCP-DENSITY-WINDOW*.cube` | Energy-window electron density | bohr⁻³ |
| `mycalc-RCP-KINETIC-ENERGY-DENSITY*.cube` | All-occupied Laplacian-form kinetic-energy density | Hartree bohr⁻³ |
| `mycalc-RCP-REGIONAL-ENERGY-DENSITY*.cube` | Energy-window regional energy density | Hartree bohr⁻³ |
| `mycalc-RCP*.cube` | Regional chemical potential | Hartree |

Filename suffixes include CP2K iteration numbers. The simple wildcard `*-RCP*.cube` matches **all four** output types; use the more specific names above to select a field.

Check the `RCP|` diagnostics for agreement between electron counts from orbital weights, the density matrix overlap trace, and real-space grid integration. Compare integrated kinetic-energy density with CP2K's kinetic energy.

**Convergence warning:** Total-energy convergence with respect to k-point sampling does not by itself imply convergence of the spatial RCP field. Compare RCP on a common grid with a consistent electron-density mask.

## Running on workstations or HPC systems

RCP does not require a particular operating-system distribution, scheduler, cluster, or MPI implementation. Start with the [single-process H₂ quick start](docs/rcp_quickstart.md). A typical patched CP2K build can then be used with your local system's MPI launch or batch scheduler if needed; follow the instructions provided for that build.

The [ISSP sham-specific notes](docs/sham_notes.md) and the optional `tools/run_cp2k_sham.slurm` template document one known HPC configuration. They **must not be treated as general CP2K requirements**.

## Reproducing validation

The numeric H2 regression references are in `examples/h2/TEST_FILES.toml`. The `tools/check_kpoint_consistency.py` script expects logs and CUBE files from several named k-point cases, laid out as specified in `examples/README.md`. **This repository does not contain those computed files.**

Original source-tree examples (the legacy five-chain large-cell restart calculation and the diamond(001) convergence benchmark) are documented in the detailed tutorial but **are not part of this distribution**. They rely on separate research datasets, analysis results, or provenance. The small, standalone examples here are the intended public entry points.

## License and attribution

The source modifications and CP2K-derived materials are distributed under **GNU GPL-2.0-or-later**, consistent with CP2K. The included [LICENSE](LICENSE) contains the GPL v2 license text; the `-or-later` permission is specified by the upstream CP2K source notices. Preserve upstream copyright and attribution notices when redistributing modifications.

This repository is an independent research extension and does not imply endorsement by the CP2K developers.

For the physical background of the RCP method, see:

> M. Fukuda *et al.*, “Regional chemical potential analysis for material surfaces,” *Journal of Chemical Physics* **164**, 084123 (2026). DOI: [10.1063/5.0288934](https://doi.org/10.1063/5.0288934).

Please also cite the CP2K software when publishing results generated with it; see [CP2K citation guidance](https://www.cp2k.org/howto:citing).

## Publishing checklist

1. Review every file and copyright/provenance note, especially before adding any new structural data.
2. Verify that the full patch applies to the stated upstream commit.
3. Keep CUBE, restart, wavefunction, logs, and binary outputs outside Git.
4. Review `git status` and `git ls-files` before pushing.
5. Update patch and documentation together when releasing a new version.
