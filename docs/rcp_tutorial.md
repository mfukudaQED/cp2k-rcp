# CP2K RCP Tutorial

> **Public-release note:** This detailed tutorial documents both the physics and selected historical validation runs from the CP2K `rcp-development` work. References to `DEV_TREE/` or historical `local-sham` datasets denote the **separate development workspace**, not material available in this public repository. Consult the [example index](../examples/README.md) for distributed inputs and the [main README](../README.md) for installable patch instructions.
>
> **Start here for routine calculations:** [CP2K RCP User Manual](rcp_user_manual.md). This longer tutorial retains implementation background, individual tests, output definitions, and practical observations.
>
> **Scope:** GPW RKS and collinear UKS, Gamma-point and general k-point sampling. **Noncollinear spin, SOC, and explicitly relativistic RCP are not supported.**

## 1. Overview

This document explains the calculation of the **regional chemical potential (RCP)** using the extension developed on the CP2K `rcp-development` branch.

Implemented and tested settings include:

- Gaussian and plane waves (GPW).
- Gamma-point calculations and general GPW k-point sampling.
- Collinear restricted and unrestricted Kohn–Sham calculations (RKS and UKS).
- Nonperiodic molecules and periodic systems.
- Integer occupations.
- Selected near-metallic examples using Fermi–Dirac smearing.

RCP is defined as

$$
\mu_R^\tau(\mathbf r)
=
\frac{\varepsilon_{\tau,\mathrm{ew}}(\mathbf r)}
     {n_{\mathrm{ew}}(\mathbf r)}.
$$

Here:

- `n_{\mathrm{ew}}(\mathbf r)` is the electron density selected by the energy window.
- `\varepsilon_{\tau,\mathrm{ew}}(\mathbf r)` is the energy-window regional energy density.
- `\mu_R^\tau(\mathbf r)` is the regional chemical potential.

The implementation can also calculate the Laplacian-form kinetic-energy density of **all occupied states**:

$$
T_e(\mathbf r)
=
-\frac14
\sum_{\mu\nu}P_{\mu\nu}
\left[
\phi_\mu(\mathbf r)\nabla^2\phi_\nu(\mathbf r)
+
\nabla^2\phi_\mu(\mathbf r)\phi_\nu(\mathbf r)
\right].
$$

This all-occupied kinetic-energy density is **not** the same as the energy-window regional energy density used in the RCP numerator.

---

## 2. Minimal input for enabling RCP

Place the RCP section at the following input hierarchy:

`DFT → PRINT → RCP`

A representative configuration is:

~~~text
&DFT
  ...
  &PRINT
    &RCP ON
      ENERGY_LOWER [eV] -3.0
      ENERGY_UPPER [eV]  0.0

      BROADENING_LOWER [eV] 0.001
      BROADENING_UPPER [eV] 0.001

      DENSITY_CUTOFF 1.0E-12

      PRINT_DENSITY_WINDOW T
      PRINT_KINETIC_ENERGY_DENSITY T
      PRINT_REGIONAL_ENERGY_DENSITY T
      PRINT_RCP T

      MPI_IO T
      STRIDE 1 1 1
    &END RCP
  &END PRINT
&END DFT
~~~

The energy window is specified **relative to the reference chemical potential**, rather than using absolute Kohn–Sham eigenvalues. For example,

~~~text
ENERGY_LOWER [eV] -3.0
ENERGY_UPPER [eV]  0.0
~~~

approximately selects states in the interval

$$
\mu-3\ {\rm eV}
<
\varepsilon
<
\mu.
$$

The exact weights are smoothly broadened at the boundaries rather than being a strict step function.

### 2.1 General k-point sampling

The k-point implementation of the CP2K 2026.2 RCP patch uses the same `&KPOINTS` input as the SCF calculation. For example, a full 2×1×1 Monkhorst–Pack mesh is specified as follows:

~~~text
&DFT
  &KPOINTS
    SCHEME MONKHORST-PACK 2 1 1
    FULL_GRID ON
    SYMMETRY OFF
    WAVEFUNCTIONS COMPLEX
  &END KPOINTS
  &PRINT
    &RCP ON
      ENERGY_LOWER [eV] -20.0
      ENERGY_UPPER [eV] 1.0
      PRINT_DENSITY_WINDOW T
      PRINT_KINETIC_ENERGY_DENSITY T
      PRINT_REGIONAL_ENERGY_DENSITY T
      PRINT_RCP T
    &END RCP
  &END PRINT
&END DFT
~~~

At each `\mathbf k`, the calculation uses complex Bloch orbitals and their eigenvalues. The energy-window density is

$$
n_{\mathrm{ew}}(\mathbf r)=
\sum_{\sigma}\sum_{\mathbf k} w_{\mathbf k}
\sum_n f_{n\mathbf k\sigma}^{\mathrm{ew}}
\left|\psi_{n\mathbf k\sigma}(\mathbf r)\right|^2,
$$

where `w_{\mathbf k}` is the normalized k-point weight supplied by CP2K. RKS includes a spin-degeneracy factor of two; UKS assigns a factor of one to each spin channel.

Temporary energy-window occupations are used to construct density matrices. The standard `kpoint_density_transform` processes Bloch phases, crystal symmetry, and real-space lattice images. The SCF occupation and k-point density-matrix state is restored around the RCP output step.

Public regression inputs are `examples/h2/H2-rcp-kpoints.inp` and `examples/h2/H2-rcp-kpoints-sym.inp` (relative to the root of this distribution). Check numerical accuracy by comparing `Tr[P_ew S]` with real-space density integration and checking the kinetic-energy integral.

---

## 3. RCP input keywords

| Keyword | Default | Description |
|---|---:|---|
| `ENERGY_LOWER` | −3.0 eV | Energy-window lower bound relative to reference chemical potential |
| `ENERGY_UPPER` | 0.0 eV | Upper energy-window bound |
| `BROADENING_LOWER` | 0.001 eV | Fermi broadening at the lower boundary |
| `BROADENING_UPPER` | 0.001 eV | Fermi broadening at the upper boundary |
| `EPS_FILTER` | 1.0E−14 | Filtering threshold for Gamma-point sparse-matrix operations |
| `DENSITY_CUTOFF` | 1.0E−12 | Cutoff for the RCP denominator `n_ew` |
| `PRINT_DENSITY_WINDOW` | T | Write the `n_ew` CUBE file |
| `PRINT_KINETIC_ENERGY_DENSITY` | T | Write the `T_e` CUBE file |
| `PRINT_REGIONAL_ENERGY_DENSITY` | T | Write the `ε_{τ,ew}` CUBE file |
| `PRINT_RCP` | T | Write the `μ_R^τ` CUBE file |
| `MPI_IO` | T | Use MPI-I/O for CUBE output |
| `STRIDE` | 1 1 1 | CUBE grid stride along X/Y/Z |

Always specify the units explicitly when entering energy-related values:

~~~text
ENERGY_LOWER [eV] -3.0
ENERGY_UPPER [eV] 0.0
BROADENING_LOWER [eV] 0.001
BROADENING_UPPER [eV] 0.001
~~~

---

## 4. Reference chemical potential

### 4.1 Calculations without smearing

For an ordinary calculation with integer occupations, the reference for the RCP energy window is

$$
\mu
=
\frac{
\varepsilon_{\rm HOMO}
+
\varepsilon_{\rm LUMO}
}{2}.
$$

In UKS, first determine across **all spin channels and k points**:

- **Global HOMO:** the largest HOMO eigenvalue among the spin channels and k points.
- **Global LUMO:** the smallest LUMO eigenvalue among those channels and k points.

The midpoint of these global levels is used. The standard output contains

~~~text
RCP| Global HOMO [a.u.]
RCP| Global LUMO [a.u.]
RCP| HOMO-LUMO midpoint [a.u.]
RCP| Reference chemical potential [a.u.]
~~~

Without smearing, **reference chemical potential = HOMO–LUMO midpoint**.

### 4.2 Fermi–Dirac smearing

When electronic smearing is active in the SCF calculation, the reference is instead the **SCF Fermi energy**, not the HOMO–LUMO midpoint:

~~~text
&SCF
  ...
  &SMEAR
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE [K] 300
  &END SMEAR
&END SCF
~~~

Thus **reference chemical potential = SCF Fermi energy**. This is the natural choice for metallic and near-metallic systems; ordinary convergence checks are still required.

---

## 5. Smooth energy-window weights

The energy window is defined by the difference between smooth Fermi functions, not by a discontinuous step function. Conceptually, each KS state has weight

$$
w_i
=
f(
\varepsilon_i-\mu;
E_{\rm upper},
\sigma_{\rm upper}
)
-
f(
\varepsilon_i-\mu;
E_{\rm lower},
\sigma_{\rm lower}
).
$$

The implementation clips the result to the interval [0,1]. A state exactly at a window boundary generally receives a fractional weight when the broadening is finite.

In a finite-molecule regression test, it may be useful to choose an `ENERGY_UPPER` slightly above 0 eV to include the HOMO almost completely. For surface and metallic systems, using

~~~text
ENERGY_UPPER [eV] 0.0
~~~

naturally selects the neighborhood below the Fermi level with a smooth boundary. Note that the window broadening is separate from the smearing temperature used to converge the parent SCF calculation.

---

## 6. Closed-shell molecular example: benzene

The standalone validated public input is

~~~text
examples/benzene/benzene.inp
~~~

and the original development input was located at `DEV_TREE/local-sham/rcp_validation/benzene/input/benzene.inp`.

Its RCP section uses:

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

    STRIDE 1 1 1
  &END RCP
&END PRINT
~~~

For a nonperiodic molecule, the `CELL` and `POISSON` periodicity settings must agree:

~~~text
&POISSON
  PERIODIC NONE
  PSOLVER WAVELET
&END POISSON

...

&CELL
  ...
  PERIODIC NONE
&END CELL
~~~

---

## 7. Open-shell collinear UKS example: C₂H₅

The present RCP implementation supports collinear UKS. The validated public input is

~~~text
examples/c2h5/c2h5.inp
~~~

and its original development location was `DEV_TREE/local-sham/rcp_validation/c2h5/input/c2h5.inp`.

For the C₂H₅ radical, use a doublet:

~~~text
&DFT
  UKS T
  MULTIPLICITY 2
  ...
&END DFT
~~~

The corresponding RCP section is:

~~~text
&PRINT
  &RCP ON
    ENERGY_LOWER [eV] -2.0
    ENERGY_UPPER [eV] 0.0

    BROADENING_LOWER [eV] 0.001
    BROADENING_UPPER [eV] 0.001

    DENSITY_CUTOFF 1.0E-12

    PRINT_DENSITY_WINDOW T
    PRINT_KINETIC_ENERGY_DENSITY T
    PRINT_REGIONAL_ENERGY_DENSITY T
    PRINT_RCP T

    STRIDE 1 1 1
  &END RCP
&END PRINT
~~~

All CUBE fields sum contributions from the alpha and beta spin channels. This implementation does **not** output separate alpha- and beta-resolved RCP fields.

---

## 8. Historical periodic, near-metallic example: five-chain system

The validated development-only input was

~~~text
DEV_TREE/local-sham/rcp_validation/5chain/input/5chain_smear_serial_full.inp
~~~

**It is not distributed with these public examples:** this large calculation depended on a separate structural dataset and SCF restart state.

In the development benchmark, the frontier gap was approximately **0.025 eV**, similar to `k_B T` at 300 K. Consequently, diagonalization with Fermi–Dirac smearing was used instead of the integer-occupation orbital-transformation (OT) method.

Representative SCF section:

~~~text
&SCF
  EPS_SCF 1.0E-7
  MAX_SCF 100
  SCF_GUESS RESTART

  ADDED_MOS 50

  &MIXING
    ALPHA 0.20
    BETA 1.0
    METHOD BROYDEN_MIXING
    NBROYDEN 12
  &END MIXING

  &SMEAR
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE [K] 300
  &END SMEAR
&END SCF
~~~

Representative RCP settings:

~~~text
&PRINT
  &RCP ON
    ENERGY_LOWER [eV] -1.0
    ENERGY_UPPER [eV] 0.0

    BROADENING_LOWER [eV] 0.01
    BROADENING_UPPER [eV] 0.01

    DENSITY_CUTOFF 1.0E-12

    PRINT_DENSITY_WINDOW F
    PRINT_KINETIC_ENERGY_DENSITY T
    PRINT_REGIONAL_ENERGY_DENSITY F
    PRINT_RCP T

    MPI_IO F
    STRIDE 1 1 1
  &END RCP
&END PRINT
~~~

Parallel MPI-I/O of the large periodic CUBE stalled during the development test. The final calculation therefore used

~~~text
MPI_IO F
~~~

to write CUBE data serially. If a large calculation stalls while writing a CUBE, trying `MPI_IO F` is a useful troubleshooting step.

---

## 9. Choosing which CUBE fields to save

### 9.1 For implementation and consistency validation

Write all four fields:

~~~text
PRINT_DENSITY_WINDOW T
PRINT_KINETIC_ENERGY_DENSITY T
PRINT_REGIONAL_ENERGY_DENSITY T
PRINT_RCP T
~~~

### 9.2 For displaying RCP on an electronic surface

If only the kinetic-density isosurface and its RCP coloring are needed, the following two CUBEs often suffice:

~~~text
PRINT_DENSITY_WINDOW F
PRINT_KINETIC_ENERGY_DENSITY T
PRINT_REGIONAL_ENERGY_DENSITY F
PRINT_RCP T
~~~

With `PRINT_RCP T`, the window electron density `n_ew` and regional energy density `ε_{τ,ew}` are still evaluated internally. Setting their output flags to `F` only avoids writing the respective CUBE files. For a quantitative comparison, consider also writing the window-density CUBE to define a robust mask.

---

## 10. Output CUBE files and physical definitions

If the `&GLOBAL` section contains

~~~text
PROJECT mycalc
~~~

the program may write the following fields. File suffixes such as `1_0` depend on CP2K iteration numbering and should not be assumed constant.

### 10.1 Energy-window electron density

~~~text
mycalc-RCP-DENSITY-WINDOW1_0.cube
~~~

Field:

$$
n_{\mathrm{ew}}(\mathbf r).
$$

Unit:

$$
{\rm bohr}^{-3}.
$$

For closed-shell RKS, the density includes the factor of two for spin degeneracy. Integrating it over the entire simulation domain gives the number of electrons selected by the energy window (up to numerical integration error).

### 10.2 Laplacian-form kinetic-energy density

~~~text
mycalc-RCP-KINETIC-ENERGY-DENSITY1_0.cube
~~~

Field:

$$
T_e(\mathbf r)
=
-\frac14
\sum_{\mu\nu}P_{\mu\nu}
\left[
\phi_\mu\nabla^2\phi_\nu
+
(\nabla^2\phi_\mu)\phi_\nu
\right].
$$

Unit:

$$
{\rm Hartree}\,{\rm bohr}^{-3}.
$$

**Important:** This all-occupied Laplacian-form kinetic-energy density is **not** the positive-definite kinetic-energy density `\tau(\mathbf r)`. Neither should it be confused with the energy-window regional energy density below.

### 10.3 Energy-window regional energy density

~~~text
mycalc-RCP-REGIONAL-ENERGY-DENSITY1_0.cube
~~~

Field:

$$
\varepsilon_{\tau,\mathrm{ew}}(\mathbf r)
=
\frac18
\sum_{\mu\nu}
P^{\mathrm{ew}}_{\mu\nu}
\left[
\phi_\mu\nabla^2\phi_\nu
+
(\nabla^2\phi_\mu)\phi_\nu
\right]
-
\frac14
\sum_{\mu\nu}
P^{\mathrm{ew}}_{\mu\nu}
\nabla\phi_\mu\cdot\nabla\phi_\nu.
$$

Unit:

$$
{\rm Hartree}\,{\rm bohr}^{-3}.
$$

The window density matrix `P^{ew}` selects the energy range specified by `ENERGY_LOWER/UPPER` and their broadenings.

### 10.4 Regional chemical potential

~~~text
mycalc-RCP1_0.cube
~~~

Field:

$$
\mu_R^\tau(\mathbf r)
=
\frac{
\varepsilon_{\tau,\mathrm{ew}}(\mathbf r)
}{
n_{\mathrm{ew}}(\mathbf r)
}.
$$

Unit:

$$
{\rm Hartree}.
$$

The CUBE title also includes

~~~text
RCP REGIONAL CHEMICAL POTENTIAL [HARTREE]
~~~

As with other CP2K CUBEs, the `1_0` portion of a filename is iteration-dependent. A broad search such as `*-RCP*.cube` will find all RCP-related fields, **not just the RCP ratio**. Prefer an exact field-specific filename pattern when selecting values for analysis.

---

## 11. DENSITY_CUTOFF

Since the RCP is a ratio,

$$
\mu_R^\tau
=
\frac{
\varepsilon_{\tau,\mathrm{ew}}
}{
n_{\mathrm{ew}}
},
$$

it is numerically vulnerable in vacuum and near orbital nodes, where `n_ew → 0`. There, small denominator values can produce very large apparent RCP values or noise.

For example, with

~~~text
DENSITY_CUTOFF 1.0E-12
~~~

the RCP is set to zero at grid points satisfying

$$
n_{\mathrm{ew}}(\mathbf r)
\le
10^{-12}.
$$

This is a **denominator mask applied when computing the RCP ratio**, not a clipping option applied only during visualization.

When evaluating global correlations, minima/maxima, or histograms, exclude extremely low-density regions from physical comparisons in addition to applying the numerical cutoff. Report the mask definition alongside quantitative results.

---

## 12. STRIDE

~~~text
STRIDE 1 1 1
~~~

writes every point of the applicable CP2K real-space grid to the CUBE.

For example,

~~~text
STRIDE 2 2 2
~~~

writes every other point in each of the X/Y/Z directions, reducing the number of written grid points approximately eightfold.

Suggested use:

- **Validation and quantitative comparison:** `STRIDE 1 1 1`.
- **Visualization of very large systems:** consider `STRIDE 2 2 2` or a larger stride if the resulting spatial resolution remains adequate.

`STRIDE` changes the sampling of the **output CUBE**, not the underlying SCF solution.

---

## 13. MPI_IO

The default is:

~~~text
MPI_IO T
~~~

Parallel CUBE output is available for many small and medium-sized systems.

In the historical five-chain periodic benchmark, however, MPI-I/O stalled during CUBE output; the successful run therefore used

~~~text
MPI_IO F
~~~

to select serial CUBE writing. If the main calculation has completed but the CUBE output remains at zero bytes or makes no progress, try `MPI_IO F`.

**This keyword affects output, not the parallelism of the SCF calculation.**

---

## 14. Diagnostics to inspect in standard output

When RCP runs, it writes lines beginning with `RCP|` to the CP2K output.

Representative labels include:

~~~text
RCP| Number of collinear spin channels:
RCP| Global HOMO [a.u.]:
RCP| Global LUMO [a.u.]:
RCP| HOMO-LUMO midpoint [a.u.]:
RCP| Reference chemical potential [a.u.]:

RCP| Relative energy window [a.u.]:
RCP| Lower/upper broadening [a.u.]:

RCP| Electron count from sum(weights):
RCP| Electron count from Tr[P_window*S]:
RCP| Electron count from grid integration:
RCP| |Tr[P*S]-grid integral|:

RCP| Integral of Laplacian-form T_e [a.u.]:
RCP| CP2K kinetic energy [a.u.]:
RCP| |Integral(T_e)-E_kin| [a.u.]:

RCP| Window Laplacian component integral [a.u.]:
RCP| Window gradient component integral [a.u.]:
RCP| Regional energy density integral [a.u.]:

RCP| Density cutoff for RCP ratio [a.u.]:
~~~

Some diagnostics are conditional on the selected output and execution path. Always check that the SCF has converged and that CP2K exited normally before interpreting the RCP diagnostics.

---

## 15. Consistency check for the window electron density

Verify that three independent electron-count calculations agree.

**Sum of Fermi/window weights:**

$$
N_{\mathrm{ew}}
=
\sum_i g_s w_i.
$$

**Density matrix and AO overlap:**

$$
N_{\mathrm{ew}}
=
{\rm Tr}
\left[
P^{\mathrm{ew}}S
\right].
$$

**Integral of real-space window density:**

$$
N_{\mathrm{ew}}
=
\int
n_{\mathrm{ew}}(\mathbf r)
d^3r.
$$

Corresponding output labels:

~~~text
Electron count from sum(weights)
Electron count from Tr[P_window*S]
Electron count from grid integration
~~~

Also check that the reported residual

~~~text
|Tr[P*S]-grid integral|
~~~

is sufficiently small for the selected basis, grid, and numerical precision. The window electron count need not be an integer and need not equal the total number of valence electrons.

---

## 16. Consistency check for kinetic-energy density

The all-occupied Laplacian-form kinetic-energy density should satisfy

$$
\int
T_e(\mathbf r)
d^3r
=
E_{\rm kin}.
$$

Inspect:

~~~text
RCP| Integral of Laplacian-form T_e [a.u.]
RCP| CP2K kinetic energy [a.u.]
RCP| |Integral(T_e)-E_kin| [a.u.]
~~~

A small final residual is a strong internal check on the real-space derivative collocation and integration. An anomalously large discrepancy warrants investigation of the numerical grid, boundary conditions, and SCF state.

---

## 17. Running CP2K with the RCP patch

### 17.1 General execution

First initialize the environment of your **patched** CP2K build:

~~~bash
source /path/to/patched-cp2k/install/cp2k_env

mpiexec -n 8 cp2k.psmp \
  -i input.inp \
  -o output.out
~~~

Use the correct number of ranks, modules, and MPI setup for your system.

### 17.2 ISSP sham

The public distribution includes a portable sham Slurm template:

~~~text
tools/run_cp2k_sham.slurm
~~~

Example from a directory containing `input.inp`:

~~~bash
sbatch --export=ALL,CP2K_ENV=/path/to/patched-cp2k/install/cp2k_env \
  /path/to/cp2k-rcp/tools/run_cp2k_sham.slurm \
  input.inp output.out
~~~

On sham, the documented working MPI launch uses **`mpiexec` instead of `srun`**, with

~~~bash
export I_MPI_COLL_EXTERNAL=no
export I_MPI_FABRICS=shm:ofi
~~~

set inside the Slurm job. A separate, pre-existing development wrapper was historically located at `DEV_TREE/local-sham/run_cp2k.slurm` and is **not** required for the public example.

---

## 18. Visualizing RCP

For surfaces, displaying RCP on an electronic interface may be more informative than viewing all values across the three-dimensional volume.

In previous validation work, the Laplacian-form kinetic-energy density isosurface

$$
T_e(\mathbf r)
=
10^{-5}
$$

was colored with the regional chemical potential

$$
\mu_R^\tau(\mathbf r).
$$

This requires two field types:

~~~text
*-RCP-KINETIC-ENERGY-DENSITY*.cube
*-RCP1_0.cube
~~~

The latter is an example iteration-specific name; adapt the suffix to your output.

For a historical comparison with OpenMX on the five-chain model, an illustrative color scale was

$$
-0.4
\le
\mu_R^\tau
\le
-0.25
\ {\rm Ha}.
$$

The original validation images and CUBEs were stored in the development workspace, `DEV_TREE/local-sham/rcp_validation/visualize/`, and are **not distributed**.

Historical output names included

~~~text
benzene_RCP.cube
benzene_Te.cube

c2h5_RCP.cube
c2h5_Te.cube

5chain_RCP.cube
5chain_Te.cube
~~~

The chosen isovalue, RCP window, mask, and color scale should be reported when publishing a visualization. The particular numerical isosurface and color range used in this historical comparison are not universal defaults.

---

## 19. CUBE whitespace compatibility

The current patched CUBE writer places explicit whitespace between successive volumetric numbers. For example,

~~~text
  0.00000E+000 -0.19538E+001 -0.20544E+001
~~~

can be parsed by viewers that cannot interpret adjoining fixed-width Fortran fields.

The public CUBE whitespace normalizer is

~~~text
tools/normalize_cube_spacing.py
~~~

To check older output:

~~~bash
python3 tools/normalize_cube_spacing.py --check old.cube
~~~

To normalize an existing CUBE:

~~~bash
python3 tools/normalize_cube_spacing.py old.cube
~~~

A directory argument recursively processes contained CUBE files. The converter **does not recompute numbers or rerun CP2K**; it preserves numeric tokens while normalizing whitespace. Check or maintain backups of important data before modifying original files.

---

## 20. Current limitations

### 20.1 k-point sampling

The implemented GPW RCP driver supports Gamma-point calculations and general k-point sampling. Complex Bloch orbitals and symmetry-reduced k points, including time-reversal/inversion-related transformations, use CP2K's standard density-matrix transformation routines.

Small regression tests cover explicit `1×1×1` and `2×1×1` meshes with and without symmetry reduction. **Independent k-point convergence is required** for larger materials, metallic systems, and other sampling meshes.

### 20.2 Collinear spin only

RCP CUBE output in this CP2K 2026.2 extension supports **RKS and collinear UKS only**. It does not automatically support post-SCF SOC spinor analysis or self-consistent noncollinear magnetism.

An isolated experimental local spinor kernel, computing charge, three magnetization components, and nonrelativistic kinetic/regional energy contributions, exists in the patched source at

~~~text
src/rcp_spinor.F
~~~

with an isolated unit test at

~~~text
tests/QS/regtest-rcp/test_spinor_kernel.f90
~~~

**This kernel is not connected to the spinor eigenstate reader, AO grid collocation, SCF machinery, or RCP CUBE output.** Accordingly, the following remain unsupported:

- Noncollinear spin states.
- Noncollinear spinor wavefunctions (as distinct from the **supported** complex Bloch orbitals of the scalar/collinear k-point code).
- Spin spirals.
- Spin–orbit coupling.

Similarly, a GPW-specific density collocator does not contain the one-center all-electron corrections needed for GAPW with ZORA/DKH. The code explicitly aborts on unsupported configurations rather than presenting incorrect all-electron or relativistic RCP values.

An internal design note for eventual noncollinear/relativistic extensions was historically kept under `docs/methods/rcp_spinor_relativistic_design.md` in the CP2K development tree. **That proposal is not a description of implemented functionality**, and its implementation is currently deferred.

**Supported spin modes:** RKS and collinear UKS.

### 20.3 Validated electronic-structure settings

The core validated settings are:

- GPW.
- GTH pseudopotentials.
- Semilocal PBE.
- Gamma-only and small k-point GPW calculations.

GAPW is **not** supported for RCP and causes an explicit runtime abort. Hybrid/HF functionals, very large general k-point meshes, and metallic smearing across arbitrary materials have not been systematically validated.

A separate development study evaluated a diamond(001) PBE slab with 300 K smearing for Gamma to `8×8` k-point meshes. Even where the total energy was nearly converged, significant differences remained in spatial RCP. Those development input files and numerical outputs are **not shipped** with this standalone distribution.

### 20.4 Large-system computational costs

After the SCF calculation, RCP performs:

- **Gamma-only path:** evaluation of `S^{-1/2}` and dense KS diagonalization.
- **k-point path:** reuse of Bloch molecular orbitals from the SCF step to construct energy-window occupations.
- Construction of the energy-window density matrix.
- Real-space orbital and orbital-derivative collocation.
- CUBE writing.

Relative to the base SCF calculation, large AO systems can require additional **memory, dense linear algebra, and I/O**.

For large calculations, consider

~~~text
MPI_IO F
~~~

if parallel CUBE I/O stalls, and

~~~text
STRIDE 2 2 2
~~~

if reduced spatial resolution is acceptable for the intended visualization.

---

## 21. Suggested RCP input blocks

These configurations illustrate starting points, **not universally converged parameters**. Vary the window, broadenings, and density cutoff for each system.

### 21.1 Molecules and insulating systems

~~~text
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

  STRIDE 1 1 1
&END RCP
~~~

### 21.2 Surface RCP visualization

~~~text
&RCP ON
  ENERGY_LOWER [eV] -3.0
  ENERGY_UPPER [eV] 0.0

  BROADENING_LOWER [eV] 0.001
  BROADENING_UPPER [eV] 0.001

  DENSITY_CUTOFF 1.0E-12

  PRINT_DENSITY_WINDOW F
  PRINT_KINETIC_ENERGY_DENSITY T
  PRINT_REGIONAL_ENERGY_DENSITY F
  PRINT_RCP T

  MPI_IO F
  STRIDE 1 1 1
&END RCP
~~~

For quantitative surface analysis, also save the window-density field to define and test the spatial mask, even if the visualization itself only uses the `T_e` isosurface and RCP values.

### 21.3 Near-metallic or metallic systems

An example SCF smearing section:

~~~text
&SMEAR
  METHOD FERMI_DIRAC
  ELECTRONIC_TEMPERATURE [K] 300
&END SMEAR
~~~

Example RCP window broadening:

~~~text
BROADENING_LOWER [eV] 0.01
BROADENING_UPPER [eV] 0.01
~~~

The SCF electronic temperature and RCP energy-window broadenings are **distinct parameters**. Optimize them for the electronic structure and numerical convergence of your system.

---

## 22. Validated example systems

### Benzene

Public input: `examples/benzene/benzene.inp`. The internal regression was stored under `DEV_TREE/local-sham/rcp_validation/benzene/`.

- Closed shell.
- Nonperiodic.
- `[-3,0]` eV energy window.

### C₂H₅ radical

Public input: `examples/c2h5/c2h5.inp`. The internal regression was stored under `DEV_TREE/local-sham/rcp_validation/c2h5/`.

- Collinear UKS doublet.
- Nonperiodic.
- `[-2,0]` eV energy window.

### Five-chain system

Historical internal data location: `DEV_TREE/local-sham/rcp_validation/5chain/`. **No standalone public input or output is included.**

- 210 carbon atoms.
- Periodic, Gamma-point calculation.
- Near-metallic electronic structure.
- Fermi–Dirac smearing at 300 K.
- `[-1,0]` eV energy window.
- `MPI_IO F` for the final CUBE output.

The separate H₂ Gamma and k-point regression inputs supplied under `examples/h2/` provide a small, practical way to validate the installation without requiring these larger development datasets.

---

## 23. Quick checklist

For a new RCP calculation, verify:

1. GPW and either Gamma-point or appropriately validated k-point sampling are used.
2. The spin formalism is RKS or **collinear** UKS.
3. The SCF calculation has converged to an appropriate tolerance.
4. The physical meaning of the chosen energy window is clear.
5. `ENERGY_LOWER < ENERGY_UPPER`.
6. `BROADENING_LOWER` and `BROADENING_UPPER` are positive.
7. `PRINT_RCP T` is enabled if RCP is needed.
8. `PRINT_KINETIC_ENERGY_DENSITY T` is enabled for electronic-interface visualization.
9. Try `MPI_IO F` if writing a large CUBE stalls.
10. Prefer `STRIDE 1 1 1` for quantitative comparisons.
11. The three window-electron-count diagnostics agree.
12. The spatial integral of `T_e` agrees with the CP2K kinetic energy.
13. Vacuum and orbital-node RCP are interpreted with the `DENSITY_CUTOFF` mask in mind.
14. Units of each CUBE field are correctly interpreted.
15. k-point convergence is established for **RCP itself**, not only for total energies.

---

## 24. Related source code

The relevant files below refer to the **patched CP2K source tree**, not to this separate distribution repository.

**RCP input definition:**

~~~text
src/input_cp2k_print_dft.F
~~~

**RCP evaluation driver:**

~~~text
src/qs_energy_window.F
~~~

**Post-SCF entry point:**

~~~text
src/qs_scf_post_gpw.F
~~~

**CUBE writer:**

~~~text
src/pw/realspace_grid_cube.F
~~~

**Additional isolated spinor-kernel unit test:**

~~~text
tests/QS/regtest-rcp/test_spinor_kernel.f90
~~~

**Historical validation dataset location (not distributed):**

~~~text
DEV_TREE/local-sham/rcp_validation/
~~~

---

## 25. References

For the theoretical background and materials-surface applications of RCP, see:

M. Fukuda *et al.*, “Regional chemical potential analysis for material surfaces,”
*Journal of Chemical Physics* **164**, 084123 (2026).
DOI: [10.1063/5.0288934](https://doi.org/10.1063/5.0288934).

For CP2K's underlying GPW and electronic-structure implementation, see the [CP2K documentation](https://manual.cp2k.org/) and [CP2K citation guidelines](https://www.cp2k.org/howto:citing).
