# CP2K RCP User Manual (GPW and Collinear Spin)

**Applies to:** CP2K 2026.2 with the RCP patch (development branch `rcp-development`).
**Updated:** 2026-10-08.

> This is a practical guide to running RCP calculations, reading the output, and verifying the results. For additional mathematical and implementation details, see the [detailed RCP tutorial](rcp_tutorial.md).
>
> **Noncollinear spin, spin–orbit coupling (SOC), and relativistic approaches such as ZORA and DKH are not currently supported by the RCP implementation.**
>
> **About the public distribution:** This repository is a standalone patch, documentation, and input-example bundle, separate from the CP2K source tree. `/path/to/cp2k-rcp` refers to a clone of this distribution, while `/path/to/patched-cp2k` refers to a separate CP2K source tree on which the patch has been applied and built. See the [main README](../README.md) for installation. Large five-chain and diamond benchmark datasets and precomputed CUBE files are not distributed.

### How to use this manual

Sections **2–3** cover the first calculation and input syntax; **4–5** explain settings and definitions; **6–8** cover CUBE files, diagnostics, and convergence; **9** identifies examples; **10** addresses problems; and **11** is a pre-run checklist.

## 1. Supported methods and prerequisites

| Calculation setting | Status | Notes |
|---|---|---|
| QUICKSTEP / GPW | Supported | Mainly validated with GTH pseudopotentials and PBE |
| Gamma point, nonperiodic molecules | Supported | Including H₂ and benzene |
| Gamma point, periodic systems | Supported | Including a historical five-chain benchmark |
| Multiple k points | Supported | Complex Bloch orbitals, symmetry reduction, and MPI k-point groups |
| Restricted Kohn–Sham (RKS) | Supported | Includes the factor of two for spin degeneracy |
| Collinear unrestricted Kohn–Sham (UKS) | Supported | Alpha and beta channels are summed |
| Fermi–Dirac smearing | Supported | Tested on a surface model at 300 K; this does not establish accuracy for all metals |
| GAPW / all-electron / ZORA / DKH | **Not supported** | Requesting RCP for these settings triggers an explicit abort |
| Noncollinear spin / SOC / spin spirals | **Not supported** | Outside the scope of the present manual |
| Hybrid functionals / Hartree–Fock | Not validated | Requires separate validation |

RCP is a **post-SCF analysis of a converged DFT electronic state**. It is not a substitute for SCF convergence or geometry optimization.

### Selecting the correct executable

After applying and compiling the patch, load your patched CP2K environment (adjust the path):

~~~bash
source /path/to/patched-cp2k/install/cp2k_env
which cp2k.psmp
~~~

An unmodified CP2K 2026.2 binary cannot be assumed to recognize this experimental RCP input. Use the **patched binary** and follow the [installation README](../README.md) and the official CP2K build instructions.

## 2. Run your first RCP calculation on sham

The following are Bash commands. Use a separate working directory so that reference inputs are not overwritten.

### 2.1 Gamma-point H₂ (RKS)

~~~bash
cd /path/to/cp2k-rcp
mkdir -p scratch/h2_gamma
cp examples/h2/H2-rcp.inp scratch/h2_gamma/
cd scratch/h2_gamma

sbatch --export=ALL,CP2K_ENV=/path/to/patched-cp2k/install/cp2k_env \
  ../../tools/run_cp2k_sham.slurm H2-rcp.inp h2_gamma.out
squeue -u "$USER"

# Inspect after the job finishes
grep 'SCF run converged' h2_gamma.out
grep 'PROGRAM ENDED AT' h2_gamma.out
grep '^ RCP|' h2_gamma.out
ls -lh *.cube
~~~

### 2.2 Multiple-k-point H₂ (RKS)

~~~bash
cd /path/to/cp2k-rcp
mkdir -p scratch/h2_kpoints
cp examples/h2/H2-rcp-kpoints.inp scratch/h2_kpoints/
cd scratch/h2_kpoints
sbatch --export=ALL,CP2K_ENV=/path/to/patched-cp2k/install/cp2k_env \
  ../../tools/run_cp2k_sham.slurm H2-rcp-kpoints.inp h2_kpoints.out
~~~

The latter input uses a full Monkhorst–Pack `2 1 1` k-point grid with complex Bloch orbitals. Both inputs are complete regression-test cases and have been run to normal completion.

**Important:** These two H₂ inputs intentionally set `PRINT_DENSITY_WINDOW F` to reduce output volume. They therefore produce **three CUBE field types**, not four. Change the setting to `T` in your copied input if you also want the window-density CUBE. With `PRINT_RCP T`, the window-density numerical integral is nevertheless reported in the standard output.

### 2.3 sham-specific MPI settings

The supplied `tools/run_cp2k_sham.slurm` uses `mpiexec` rather than `srun` to launch MPI. Its essential execution commands are:

~~~bash
export I_MPI_COLL_EXTERNAL=no
export I_MPI_FABRICS=shm:ofi
export OMP_NUM_THREADS=1
mpiexec -n "$SLURM_NTASKS" cp2k.psmp -i input.inp -o output.out
~~~

Run these commands **inside a Slurm allocation**. Adjust the number of MPI ranks, OpenMP threads, and wall time through the Slurm options. The excerpt illustrates a simple two-rank, one-thread-per-rank setup.

## 3. Writing a CP2K input

### 3.1 Enabling RCP

Insert `&RCP ON` under `&FORCE_EVAL / &DFT / &PRINT` in a normal CP2K input:

~~~text
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    ! Add the normal BASIS_SET_FILE_NAME, POTENTIAL_FILE_NAME,
    ! MGRID, SCF, XC, and other DFT settings here.
    &PRINT
      &RCP ON
        ENERGY_LOWER [eV] -3.0
        ENERGY_UPPER [eV] 0.0
        BROADENING_LOWER [eV] 0.05
        BROADENING_UPPER [eV] 0.05
        DENSITY_CUTOFF 1.0E-7

        PRINT_DENSITY_WINDOW T
        PRINT_KINETIC_ENERGY_DENSITY T
        PRINT_REGIONAL_ENERGY_DENSITY T
        PRINT_RCP T

        MPI_IO F
        STRIDE 1 1 1
      &END RCP
    &END PRINT
  &END DFT
&END FORCE_EVAL
~~~

**This fragment only shows where to place the RCP section; it is not a complete runnable input.** Complete examples are available in `examples/h2/`, `examples/benzene/`, and `examples/c2h5/`.

### 3.2 RKS and UKS

Example of a collinear UKS doublet:

~~~text
&DFT
  UKS T
  MULTIPLICITY 2
  ...
&END DFT
~~~

The validated open-shell molecular example is `examples/c2h5/c2h5.inp`. The k-point UKS code-path regression input is `examples/h2/H2-rcp-kpoints-uks.inp` (an artificial UKS singlet).

**All four RCP CUBE field types contain sums over alpha and beta spin channels.** They are not spin-resolved RCP output.

### 3.3 Gamma point versus multiple k points

Omit `&KPOINTS` for a Gamma-point-only calculation. For example, to request a `4×4×1` grid:

~~~text
&KPOINTS
  SCHEME MONKHORST-PACK 4 4 1
  FULL_GRID OFF
  SYMMETRY ON
  WAVEFUNCTIONS COMPLEX
&END KPOINTS
~~~

- `FULL_GRID OFF` with `SYMMETRY ON`: use crystal symmetry to reduce the k-point set.
- `FULL_GRID ON` with `SYMMETRY OFF`: use the entire grid, for example as an independent reference.
- `PARALLEL_GROUP_SIZE 1`: a regression case exercises the k-point MPI-group path.
- An explicitly specified `1 1 1` grid and an omitted `&KPOINTS` section follow distinct internal execution paths. Their RCP values have been checked for numerical consistency, but **bitwise agreement is not guaranteed**.

### 3.4 Simulation cell and pseudopotentials

For isolated molecules, set `&CELL / PERIODIC NONE` consistently with `&POISSON / PERIODIC NONE`. For slabs, an illustrative approach is a three-dimensionally periodic supercell with sufficient vacuum and a `k_x k_y 1` sampling mesh. Select the actual Poisson solver, vacuum thickness, and dipole correction according to ordinary CP2K electrostatics requirements.

With GPW and GTH pseudopotentials, the density used for RCP is the **valence-electron density**, not the all-electron density of GAPW.

## 4. Input keywords and suggested settings

The following defaults have been checked against `src/input_cp2k_print_dft.F` in the patched CP2K source. Place every keyword under `&DFT / &PRINT / &RCP`.

| Keyword | Default | Meaning |
|---|---:|---|
| `ENERGY_LOWER` | −3.0 eV | Lower energy-window boundary relative to reference chemical potential |
| `ENERGY_UPPER` | 0.0 eV | Upper energy-window boundary relative to the same reference |
| `BROADENING_LOWER` | 0.001 eV | Fermi-function broadening at the lower boundary |
| `BROADENING_UPPER` | 0.001 eV | Fermi-function broadening at the upper boundary |
| `EPS_FILTER` | 1.0×10⁻¹⁴ | Filtering threshold used in Gamma-point sparse-matrix operations |
| `DENSITY_CUTOFF` | 1.0×10⁻¹² a.u. | Threshold for the window-density denominator of the RCP ratio |
| `PRINT_DENSITY_WINDOW` | T | Write the energy-window electron-density CUBE |
| `PRINT_KINETIC_ENERGY_DENSITY` | T | Write the all-occupied Laplacian-form kinetic-energy density |
| `PRINT_REGIONAL_ENERGY_DENSITY` | T | Write the energy-window regional energy density |
| `PRINT_RCP` | T | Write the RCP CUBE |
| `MPI_IO` | T | Use MPI-I/O to write CUBE files |
| `STRIDE` | 1 1 1 | Sampling stride in the X, Y, Z directions |

Set the energy window explicitly rather than relying on defaults; **always specify the energy unit, preferably `[eV]`**.

### Practical example settings

| Purpose | Example window | Example broadening | Density cutoff | Outputs |
|---|---|---|---|---|
| Narrow window in a molecule/insulator | −3 to 0 eV | 0.001 eV | 1×10⁻¹² | All four fields; `STRIDE 1 1 1` |
| Surface RCP analysis | −3 to 0 eV | 0.05 eV | 1×10⁻⁷ | Density, RCP, regional energy; kinetic density as required |
| Near-metallic exploration | −1 to 0 eV | 0.01 eV | System-dependent | Start by saving both window density and RCP |
| Reduce CUBE storage | Application-dependent | As appropriate | As appropriate | `STRIDE 2 2 2` and `F` for unwanted fields |

**These are illustrative values, not universal recommendations.** Converge the RCP field with respect to window width, broadening, and masking threshold.

- `ENERGY_LOWER` must be strictly smaller than `ENERGY_UPPER`.
- Both `BROADENING_LOWER` and `BROADENING_UPPER` must be positive.
- With `PRINT_RCP T`, the numerator and denominator of the RCP ratio are evaluated internally even if their individual CUBE output flags are `F`.
- Increasing `DENSITY_CUTOFF` masks more numerically unstable near-vacuum and orbital-node regions by setting RCP to zero. **It is not merely a plotting threshold.**
- `EPS_FILTER` primarily controls filtering during the Gamma-point density-matrix construction and orthogonalization; it does not set the window occupations in the k-point implementation.
- `MPI_IO F` changes the **CUBE output method**; it does not disable MPI execution of the SCF calculation.
- `STRIDE` subsamples the written CUBE fields, but does not recompute the SCF or the internal RCP field on a coarser grid.

### Before choosing the window

Distinguish **which electronic states you wish to select** from **where you wish to display the resulting RCP**. To relate dangling-bond distributions to adsorption sites, inspect the window electron density first to ensure the relevant orbitals are included, then visualize the RCP and kinetic-energy density.

## 5. Definition of RCP and window occupations

### 5.1 Reference energy

Window boundaries are measured **relative to a reference chemical potential**, not from an absolute Kohn–Sham eigenvalue origin:

```math
E_{\mathrm{lower}} \le \varepsilon_{n\mathbf k\sigma}-\mu \le E_{\mathrm{upper}}.
```

- **No SCF smearing:** form `μ=(ε_HOMO+ε_LUMO)/2` from the global HOMO and LUMO across all spin channels and k points.
- **Fermi–Dirac SCF smearing:** use the Fermi energy determined by CP2K's SCF calculation as `μ`.

Thus `ENERGY_LOWER [eV] -3.0` and `ENERGY_UPPER [eV] 0.0` select orbital components approximately between −3 eV and 0 eV relative to the reference.

### 5.2 Smooth energy-window weights

The window selects each KS eigenstate by subtracting two smoothed Fermi functions and restricting the resulting weight to the interval [0, 1]:

```math
F(x;E,\delta)=\frac{1}{1+\exp[(x-E)/\delta]}.
```

```math
\begin{aligned}
D_{n\mathbf{k}\sigma}
&=F(\varepsilon_{n\mathbf{k}\sigma}-\mu;E_{\mathrm{upper}},\delta_{\mathrm{upper}})\\
&\quad-F(\varepsilon_{n\mathbf{k}\sigma}-\mu;E_{\mathrm{lower}},\delta_{\mathrm{lower}}),\\
W_{n\mathbf{k}\sigma}
&=\min\left\{1,\max\left\{0,D_{n\mathbf{k}\sigma}\right\}\right\}.
\end{aligned}
```

The clipping operation restricts values numerically to [0,1]. **SCF electronic-temperature smearing and RCP's `BROADENING_LOWER/UPPER` are independent.** States near the window edges receive fractional window weights; a noninteger window electron count is therefore not inherently an error.

### 5.3 Real-space quantities

The spin-summed energy-window electron density is

```math
n_{\mathrm{EW}}(\mathbf r)=
\sum_{\sigma,n,\mathbf k}
w_{\mathbf k}\,g_{\sigma}W_{n\mathbf k\sigma}\,
|\psi_{n\mathbf k\sigma}(\mathbf r)|^2.
```

For RKS, `g=2` accounts for spin degeneracy; each spin channel in UKS has `g=1`.

Given the energy-window density matrix `P_EW`, define the Laplacian-form and gradient-form kinetic-energy densities:

```math
t_L^{\rm EW}(\mathbf r)
=-\frac14\sum_{\mu\nu}P^{\rm EW}_{\mu\nu}
\bigl[\phi_\mu\nabla^2\phi_\nu+(\nabla^2\phi_\mu)\phi_\nu\bigr],
```

```math
t_G^{\rm EW}(\mathbf r)
=\frac12\sum_{\mu\nu}P^{\rm EW}_{\mu\nu}
\,\nabla\phi_\mu\cdot\nabla\phi_\nu.
```

The implemented regional energy density and RCP are

```math
\varepsilon_{\tau,\rm EW}(\mathbf r)
=-\frac12\left[t_L^{\rm EW}(\mathbf r)+t_G^{\rm EW}(\mathbf r)\right],
\qquad
\mu_R^\tau(\mathbf r)
=\frac{\varepsilon_{\tau,\rm EW}(\mathbf r)}{n_{\rm EW}(\mathbf r)}.
```

**Grid points with `n_EW ≤ DENSITY_CUTOFF` have their RCP set to zero.**

The output `T_e` from `PRINT_KINETIC_ENERGY_DENSITY` is the Laplacian-form kinetic-energy density for **all occupied states**, not the window-selected `t_L^EW`. Do not confuse this field with the numerator of the RCP ratio.

## 6. Reading results and CUBE files

For `&GLOBAL / PROJECT mycalc`, representative filenames are shown below. Suffixes such as `1_0` depend on CP2K iteration numbering; **do not hard-code them into analysis scripts**.

| Output keyword | Representative CUBE file | Field | Unit |
|---|---|---|---|
| `PRINT_DENSITY_WINDOW T` | `mycalc-RCP-DENSITY-WINDOW1_0.cube` | Window electron density `n_EW(r)` | bohr⁻³ |
| `PRINT_KINETIC_ENERGY_DENSITY T` | `mycalc-RCP-KINETIC-ENERGY-DENSITY1_0.cube` | All-occupied Laplacian-form `T_e(r)` | Ha bohr⁻³ |
| `PRINT_REGIONAL_ENERGY_DENSITY T` | `mycalc-RCP-REGIONAL-ENERGY-DENSITY1_0.cube` | `ε_{τ,EW}(r)` | Ha bohr⁻³ |
| `PRINT_RCP T` | `mycalc-RCP1_0.cube` | `μ_R^τ(r)` | Ha |

The RCP field is in **Hartree**. Multiply by approximately **27.211386** to obtain eV. CUBE header coordinates and grid vectors conventionally use bohr; always confirm the header conventions of the generated file.

The four fields from the same SCF calculation share a real-space grid. Combine only fields corresponding to the **same calculation, `STRIDE`, and iteration index**.

### 6.1 Common visualization and analysis approaches

- **Map the spatial RCP:** Load the `*-RCP1_0.cube` field into a suitable volumetric-data viewer.
- **Color an electronic surface by RCP:** Create an isosurface of `T_e(r)` and color it using `μ_R^τ(r)`. Report the isovalue and color scale.
- **Assess the reliability of the ratio:** Export `n_EW(r)` and exclude vacuum and orbital-node regions with very small window densities.
- **Compare local RCP with electron density:** Display fields using the same coordinate system and avoid drawing conclusions solely from global minima or maxima.

### 6.2 CUBE whitespace and viewer compatibility

The current RCP code writes explicit whitespace between successive numerical CUBE values. If older files cannot be read by some viewers, the existing values can be reformatted **without rerunning CP2K**:

~~~bash
cd /path/to/cp2k-rcp
python3 tools/normalize_cube_spacing.py --check /path/to/cube_file.cube
python3 tools/normalize_cube_spacing.py --dry-run /path/to/cube_directory
python3 tools/normalize_cube_spacing.py /path/to/cube_directory
~~~

Passing a directory recursively processes CUBE files beneath it. Confirm that backups of valuable original CUBE data exist before using the converter.

## 7. Validating the standard output

Inspect the following lines after execution:

~~~bash
grep 'SCF run converged' output.out
grep 'PROGRAM ENDED AT' output.out
grep '^ RCP|' output.out
~~~

**Do not use an RCP field from an unconverged SCF calculation as the basis for a physical conclusion.** Depending on the enabled fields, CP2K reports the following diagnostics. Each label is preceded by the literal output prefix `RCP|`:

| Diagnostic label | Interpretation |
|---|---|
| `Global HOMO`, `Global LUMO` | Band edges used to construct the energy reference |
| `HOMO-LUMO midpoint` | Energy reference when SCF smearing is absent |
| `Reference chemical potential` | Reference `μ` actually used for the window (Ha) |
| `Relative energy window` | The input relative bounds, printed in Ha |
| `Electron count from sum(weights)` | Sum of window occupations over KS states and k-point weights |
| `Electron count from Tr[P_window*S]` | Electron count from the window density and overlap matrices |
| `Electron count from grid integration` | Real-space integral of the window density |
| `Tr[P*S]-grid integral` (absolute difference) | Consistency of transformation, collocation, and integration |
| `Integral of Laplacian-form T_e` | Real-space integral of the all-occupied Laplacian kinetic-energy density (Ha) |
| `CP2K kinetic energy` | Kinetic energy from the CP2K electronic-structure calculation (Ha) |
| `Integral(T_e)-E_kin` (absolute difference) | Independent kinetic-energy consistency diagnostic |
| `Window Laplacian component integral` | Integral of the energy-window Laplacian component |
| `Window gradient component integral` | Integral of the energy-window gradient component |
| `Regional energy density integral` | Volume integral of the window regional energy density |

The two absolute-difference labels include vertical bars in the actual CP2K output:

~~~text
RCP| |Tr[P*S]-grid integral|:
RCP| |Integral(T_e)-E_kin| [a.u.]:
~~~

### 7.1 Three independent electron-count checks

```math
N_{\rm EW}^{\rm weights}
\stackrel{?}{=}
\mathrm{Tr}[P_{\rm EW}S]
\stackrel{?}{=}
\int n_{\rm EW}(\mathbf r)\,d^3r .
```

**Check agreement between all three evaluations.** The number of electrons in the selected energy window need not equal the total number of valence electrons.

For a validated H₂ example whose window included essentially all occupied valence states:

~~~text
RCP| Electron count from grid integration: 2.000000002782
RCP| Regional energy density integral [a.u.]: -1.134090736496
~~~

These values **depend on the system, structure, basis, k-point sampling, window, and computational settings** and are not general target values.

### 7.2 Independent kinetic-energy check

```math
\int T_e(\mathbf r)\,d^3r \stackrel{?}{=} E_{\rm kinetic}^{\rm CP2K}.
```

For a complete integral with suitable boundary conditions, the window Laplacian- and gradient-form kinetic-energy integrals should also agree. Investigate the SCF state, integration grid, boundary conditions, and k-point transformation if the discrepancy is unusually large.

### 7.3 Missing diagnostics

Check that `&RCP ON` appears in the correct `&DFT / &PRINT` section, the output verbosity (`PRINT_LEVEL`), the SCF completion, and the version of the executable. The presence of `RCP|` lines **does not by itself demonstrate SCF convergence**.

## 8. k-point convergence and surface analysis

### 8.1 Conditions to keep fixed when changing the mesh

- **Identical atomic geometry and unit cell**; do not relax the structure separately for each mesh.
- Same functional, basis, pseudopotential, `CUTOFF`, and `REL_CUTOFF`.
- Same `SCF EPS_SCF`, SCF smearing, energy-window bounds, and window broadening.
- Same real-space grid, `STRIDE`, and `DENSITY_CUTOFF`.
- Same analysis plane and electron-window density mask.

Track **more than the total energy**: compare window electron count, spatial window-density distributions, RCP cross sections/isosurfaces, and RMS differences. Comparison between k meshes requires consistent real-space coordinates.

### 8.2 Quantitative comparisons excluding vacuum and nodes

The ratio `μ_R^τ = ε_{τ,EW}/n_EW` becomes unstable when `n_EW` approaches zero. Construct a reference mask `Ω` using the window density from the reference calculation and use **the same set of grid points** for all meshes.

```math
\Delta_{\rm RMS}(\mathcal{K},\mathcal{K}_{\mathrm{ref}})
=
\sqrt{\frac{1}{|\Omega|}
\sum_{\mathbf r_i\in\Omega}
\left[
\mu_R^{\tau,\mathcal{K}}(\mathbf r_i)
-\mu_R^{\tau,\mathcal{K}_{\mathrm{ref}}}(\mathbf r_i)
\right]^2}.
```

Test mask thresholds such as 0.1% and 5% of the maximum window density. The densest calculated k mesh is only a **provisional reference**; a zero RMS difference to itself is not evidence that the physical RCP field has converged.

### 8.3 Previous diamond(001) benchmark

An internal development study compared surface-slab calculations from Gamma-only to `8×8` sampling. **The input files, CUBE fields, analysis graphics, and original structural data are not distributed** with this public release. The following numbers illustrate convergence issues, not a result independently reproducible from this repository alone.

The total-energy difference between `6×6` and `8×8` was approximately **1.84 meV per cell**, whereas the surface-RCP RMS difference remained approximately **3.89 eV**. **RCP k-point convergence was not demonstrated through 8×8.** The fixed-geometry model contained C₁₀H₄ with approximate bottom hydrogen termination. Quantitative physical predictions would also require reconsidering slab thickness and atomic relaxation.

## 9. Included examples and validation resources

The following paths are relative to the public distribution root `/path/to/cp2k-rcp`.

| Purpose | Input or validation file |
|---|---|
| H₂, isolated, Gamma-point RKS | `examples/h2/H2-rcp.inp` |
| H₂, periodic Gamma-point RKS | `examples/h2/H2-rcp-gamma-periodic.inp` |
| H₂, explicit `1×1×1` k mesh | `examples/h2/H2-rcp-kpoints-1x1x1.inp` |
| H₂, `2×1×1` full k grid | `examples/h2/H2-rcp-kpoints.inp` |
| H₂, `2×1×1` symmetry-reduced k grid | `examples/h2/H2-rcp-kpoints-sym.inp` |
| H₂, k-point UKS singlet | `examples/h2/H2-rcp-kpoints-uks.inp` |
| H₂, k-point MPI groups | `examples/h2/H2-rcp-kpoints-parallel.inp` |
| Benzene molecule, RKS | `examples/benzene/benzene.inp` |
| C₂H₅ radical, UKS doublet | `examples/c2h5/c2h5.inp` |

### 9.1 Checking the k-point regression results

If you have run the H₂ k-point inputs and assembled the required `*.out` and `*.cube` outputs in one `regression-results/` directory, check their SCF diagnostics, real-space integrations, and symmetry consistency using:

~~~bash
cd /path/to/cp2k-rcp
python3 tools/check_kpoint_consistency.py regression-results
~~~

**Precomputed output files are not shipped.** The reference values of the existing CP2K regression cases are recorded in `examples/h2/TEST_FILES.toml`. See [examples/README.md](../examples/README.md) for the required filenames and cases.

### 9.2 Preserving your calculation metadata

Record at least the following with your own results:

- CP2K input `*.inp`, standard output `*.out`, and SCF convergence status.
- CP2K executable/build revision and patch revision.
- Basis sets, pseudopotentials, functional, lattice, and atomic coordinates.
- k-point sampling, Fermi smearing, window bounds, and window broadening.
- CUBE files, `DENSITY_CUTOFF`, `STRIDE`, and masks used in plotting/comparisons.

## 10. Troubleshooting

| Symptom | Check and remedy |
|---|---|
| `RCP` input section unrecognized | Use a `cp2k.psmp` binary compiled with the patch, not unmodified CP2K. |
| `SCF run NOT converged` | Fix the SCF first. Review `EPS_SCF`, `MAX_SCF`, `MIXING`, `ADDED_MOS`, smearing, and ordinary DFT settings. |
| `RCP ENERGY_LOWER must be smaller ...` | Require `ENERGY_LOWER < ENERGY_UPPER`. |
| `RCP broadenings must be positive` | Set both window broadenings to positive values. |
| `RCP ... needs occupied and unoccupied bands` | The k-point path needs occupied and unoccupied bands; check `ADDED_MOS` and the available eigenstates. |
| Huge RCP values in vacuum | Check window density and `DENSITY_CUTOFF`; exclude very-low-density regions from quantitative comparison. |
| CUBE files are too large | Try `STRIDE 2 2 2`, disable unnecessary `PRINT_*` flags, or restrict the analysis area. |
| CUBE writing hangs or a viewer cannot parse a file | Try `MPI_IO F`. Check legacy files using `normalize_cube_spacing.py`. |
| Three window electron counts disagree | Check SCF, k-point settings, grid integration, and density-matrix transformations. |
| `GAPW ... one-center ... not implemented` | The current RCP output is GPW-only; all-electron GAPW RCP cannot be produced. |
| `RCP ... relativistic kinetic operator ... not implemented` | Explicit relativistic kinetic operators are not supported. |
| k-point results differ substantially from Gamma-only | Check geometry, periodic boundaries, meshes, electron counts, and convergence of the RCP field itself. |
| `srun` gives an MPI error | Use `mpiexec` and the designated `I_MPI_*` environment settings on sham. |

Do not conceal a large discrepancy merely by increasing `DENSITY_CUTOFF`. Examine both the window density and regional energy density to determine why the RCP ratio is sensitive.

## 11. Pre-run checklist

1. [ ] Use a **patched CP2K 2026.2 executable**.
2. [ ] Use GPW with RKS or collinear UKS.
3. [ ] Choose a physically sound cell, boundary conditions, Poisson solver, basis, and pseudopotentials.
4. [ ] Specify Gamma-only or k-point sampling and plan k-point convergence tests.
5. [ ] Ensure the SCF converges and enough unoccupied bands are available if required.
6. [ ] Understand the reference `μ`, `ENERGY_LOWER/UPPER`, and broadening.
7. [ ] Consider `PRINT_DENSITY_WINDOW T` to validate the spatial mask used with RCP.
8. [ ] Remember that `PRINT_KINETIC_ENERGY_DENSITY T` gives the **all-occupied** `T_e`.
9. [ ] Examine the `RCP|` electron-count and kinetic-energy checks.
10. [ ] Record RCP units (Ha), cutoff, `STRIDE`, isovalue, and color scale.
11. [ ] Demonstrate convergence of the RCP distribution itself, not just total energy.
12. [ ] On sham, use `mpiexec` with the required MPI environment variables.

---

**Further reading and source references**

- [Detailed RCP tutorial](rcp_tutorial.md): theory, molecular/periodic examples, and visualization.
- [Public example inputs](../examples/README.md) and `examples/h2/` regression reference values.
- Diamond(001) benchmark data are held separately and are **not included in this distribution**.
- [Installation and patch instructions](../README.md).
- `src/input_cp2k_print_dft.F` in the patched CP2K source: `&PRINT / &RCP` keyword definitions.
- `src/qs_energy_window.F` in the patched CP2K source: RCP implementation and Gamma/k-point density-matrix paths.
