# RCP with CP2K: A Beginner's Quick Start

**Audience:** Readers who are new to density-functional theory (DFT), scientific-computing clusters, or CP2K. This guide uses a small hydrogen molecule and a terminal on an ordinary Linux/macOS workstation or other computer where CP2K is available. **No Slurm scheduler, MPI knowledge, or university computing cluster is required.**

**Estimated task:** Run an already-prepared example and check its outputs. Installing or compiling CP2K can take additional time and may require help from a system administrator.

## 1. What is RCP?

Atoms and molecules contain electrons. Different electronic states contribute to how a molecule interacts with its surroundings. **Regional chemical potential (RCP)** is a way to study *where* selected electrons may make one region of a system different from another.

This CP2K extension calculates RCP on a three-dimensional grid after a normal electronic-structure calculation has finished. It can write grid files that show the selected electron density, a regional energy density, and the resulting RCP. You do **not** need to understand the mathematics to run this tutorial.

Three words you will encounter:

- **Input file (`.inp`):** plain text describing the atoms and the calculation.
- **Output log (`.out`):** text explaining what CP2K did and whether it converged (i.e., settled on a sufficiently stable answer).
- **CUBE file (`.cube`):** a three-dimensional grid of numbers that a visualization program can display.

RCP is a research method. An attractive color plot is **not, by itself, proof of adsorption preference or chemical reactivity**.

## 2. What you need

1. A computer with a terminal and a **working executable of CP2K 2026.2 built with the RCP patch**. The normal, unmodified CP2K release will not recognize the `&RCP` input section.
2. CP2K's standard `BASIS_MOLOPT` and `GTH_POTENTIALS` data files, installed where CP2K can find them.
3. The files in this public repository, especially `examples/h2/H2-rcp.inp`.
4. A text editor. A scientific CUBE viewer is optional.

**Already have a patched CP2K installation?** Go directly to Step 3.

**Do not have one?** Follow the [main README](../README.md) to download CP2K 2026.2, apply the patch, and build CP2K using its official installation instructions. The precise build commands depend on your operating system and compiler. On a managed cluster, ask the administrator how to load your patched CP2K build. Do not assume an existing `cp2k.psmp` command is patched.

The examples below use Bash syntax, available on Linux and macOS and in Windows Subsystem for Linux (WSL). On Windows without WSL, adapt the shell commands to your environment.

## 3. Prepare a clean working directory

Open a terminal. First, go to the folder into which you cloned or extracted **this** `cp2k-rcp` repository:

~~~bash
cd /path/to/cp2k-rcp
mkdir -p scratch/my_first_rcp
cp examples/h2/H2-rcp.inp scratch/my_first_rcp/
cd scratch/my_first_rcp
~~~

Replace `/path/to/cp2k-rcp` with the actual directory on **your** computer.

This keeps the original input unchanged. The `scratch/` directory is ignored by Git and is not uploaded to GitHub.

## 4. Confirm CP2K is available

If the executable is on your `PATH`, type:

~~~bash
command -v cp2k.psmp
~~~

You should see a filesystem path. If nothing is printed, ask your administrator how to load the patched binary, or use its absolute path.

A build may provide `cp2k.ssmp` instead of `cp2k.psmp`, depending on the compilation configuration. Either works for a single-process introductory run **if the binary contains this RCP extension**.

If CP2K cannot locate `BASIS_MOLOPT` or `GTH_POTENTIALS`, the input is not necessarily wrong: configure CP2K's data-file search path or replace these names in a *copy* of the input with valid paths to the installed CP2K data files. See the [CP2K manual](https://manual.cp2k.org/) for your build's data-file setup.

## 5. Run the H₂ example

From `scratch/my_first_rcp/`:

~~~bash
cp2k.psmp -i H2-rcp.inp -o H2-rcp.out
~~~

For a build providing only `cp2k.ssmp`, use `cp2k.ssmp` instead. You do not need `mpiexec` or a batch scheduler for this first calculation.

CP2K calculates the electronic structure of a hydrogen molecule and then evaluates the RCP grid. The process may take a little time. Output files are created in the **current directory**.

When it finishes, check:

~~~bash
grep 'SCF run converged' H2-rcp.out
grep 'PROGRAM ENDED AT' H2-rcp.out
grep 'RCP|' H2-rcp.out
ls -lh *.cube
~~~

**What should you see?**

- `SCF run converged` means that the *self-consistent field* (SCF) calculation—the repeated steps used to find the electrons' distribution—reached its stopping criterion.
- `PROGRAM ENDED AT` means CP2K finished normally.
- Lines beginning with `RCP|` provide diagnostic values from the RCP analysis.
- The `*.cube` files contain three-dimensional data.

A previously validated execution of this input gave a window electron count close to **2**, as expected when essentially all two H₂ electrons are selected. Small numerical deviations are normal. Exact energies depend on the input, executable, and numerical settings.

## 6. Understand the resulting files

The provided H₂ input uses project name `H2-rcp-regtest`. Its normal default output includes **three CUBE field types**:

| Example filename pattern | Meaning | Unit |
|---|---|---|
| `H2-rcp-regtest-RCP1_0.cube` | Regional chemical potential, RCP | Hartree |
| `H2-rcp-regtest-RCP-KINETIC-ENERGY-DENSITY1_0.cube` | Kinetic-energy density for all occupied states | Hartree per bohr³ |
| `H2-rcp-regtest-RCP-REGIONAL-ENERGY-DENSITY1_0.cube` | Energy-window regional energy density | Hartree per bohr³ |

The numerical suffix may differ. These are **data**, not image files. A CUBE-capable visualization tool can display slices, surfaces, or color maps.

The fourth optional field is the **energy-window electron density**, enabled by changing

~~~text
PRINT_DENSITY_WINDOW F
~~~

to

~~~text
PRINT_DENSITY_WINDOW T
~~~

in the **copied** input file and running CP2K again. Its name contains `RCP-DENSITY-WINDOW` and its unit is bohr⁻³.

**Caution:** The provided H₂ input uses `STRIDE 4 4 4`, a coarsely sampled grid chosen to make regression tests quick and compact. It is useful to verify that the program works, **not** for high-quality images or quantitatively converged spatial comparisons.

## 7. Change just one RCP setting

In the input file, find the section

~~~text
&PRINT
  &RCP ON
    ...
  &END RCP
&END PRINT
~~~

The following settings have simple meanings:

| Input keyword | Meaning |
|---|---|
| `ENERGY_LOWER` and `ENERGY_UPPER` | Which electronic energy range to select, measured in eV relative to a reference energy |
| `BROADENING_LOWER` and `BROADENING_UPPER` | How gradually the selection changes near each end |
| `PRINT_RCP T` | Save the RCP CUBE |
| `PRINT_DENSITY_WINDOW T` | Also save the selected electron-density CUBE |
| `STRIDE 4 4 4` | Write every fourth grid point in each direction; `1 1 1` writes every point |
| `DENSITY_CUTOFF` | Avoid an unstable division in regions with almost no selected electron density |
| `MPI_IO F` | Choose serial writing for the output CUBEs, even when the main calculation is parallel |

For a first change, set `PRINT_DENSITY_WINDOW T`. If you also change `STRIDE 4 4 4` to `STRIDE 1 1 1`, the generated CUBEs become larger and spatially denser. Keep the original values until the example has run successfully.

**Do not copy a different energy window blindly into a scientific project.** Energy ranges and numerical settings that make sense for H₂ may be inappropriate for an extended surface or metal.

## 8. Common problems

| Problem | What to check |
|---|---|
| `cp2k.psmp: command not found` | The CP2K executable is not on your `PATH`. Load the appropriate environment or use an absolute executable path. |
| CP2K says `RCP` is an unknown input section | You are probably using **unpatched** CP2K. Recheck the build. |
| CP2K cannot find `BASIS_MOLOPT` or `GTH_POTENTIALS` | Point the input or data-file configuration to your CP2K data directory. |
| There is no `SCF run converged` line | Do not interpret the RCP result. Inspect the output log for SCF errors. |
| No CUBE files appear | Check that CP2K finished, `&RCP ON` is present, and the CUBE output flags are `T`. |
| The CUBE is difficult to visualize | Verify you chose the RCP field (not kinetic energy), and check the plotting range and the density cutoff. |
| The CUBE looks blocky | This H₂ test uses `STRIDE 4 4 4`. Try `STRIDE 1 1 1` after the first successful run. |

## 9. Where to go next

For a real research calculation, the input must describe your **own physical system** correctly, and you must check SCF, basis/grid, k-point, energy-window, and RCP convergence. A calculation that runs successfully is not necessarily physically reliable.

- [RCP User Manual](rcp_user_manual.md): routine calculations, settings, units, diagnostics, and convergence.
- [Detailed RCP Tutorial](rcp_tutorial.md): mathematical definitions and worked research examples.
- [Example inputs](../examples/README.md): H₂, benzene, and a C₂H₅ radical.
- [Main README](../README.md): patch installation and scope.

**Unsupported:** GAPW/all-electron RCP, explicit relativistic RCP, spin–orbit coupling, and noncollinear spin remain outside the validated features.
