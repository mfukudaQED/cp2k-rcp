# Site-specific example: ISSP sham (optional)

**This is an optional platform note, not a prerequisite for RCP.** The main [Quick Start](rcp_quickstart.md), [User Manual](rcp_user_manual.md), and [Detailed Tutorial](rcp_tutorial.md) use a generic CP2K environment and can be followed without a batch scheduler or an ISSP account.

The example below is only for the ISSP sham HPC environment and the Intel MPI setup on which the RCP development tests were run.

## Environment and MPI launch

On that particular system, the working launch configuration was:

~~~bash
export I_MPI_COLL_EXTERNAL=no
export I_MPI_FABRICS=shm:ofi
~~~

MPI execution is launched through `mpiexec` rather than `srun` because of site-specific compatibility issues. **Do not copy these settings to another computer unless its administrator specifically recommends them.**

The distribution includes `tools/run_cp2k_sham.slurm`. This is an optional Slurm template, not a general installation requirement. It expects:

- An installed patched CP2K executable named `cp2k.psmp`;
- the patched build's environment initialization file, passed as `CP2K_ENV`;
- Slurm, Intel MPI, and the specific module settings available on sham.

For example, from a directory containing `H2-rcp.inp`:

~~~bash
sbatch --export=ALL,CP2K_ENV=/path/to/patched-cp2k/install/cp2k_env \
  /path/to/cp2k-rcp/tools/run_cp2k_sham.slurm \
  H2-rcp.inp H2-rcp.out
~~~

Choose appropriate resources and paths for the actual job. On **other clusters**, use the scheduler, resource allocation, MPI launch, and environment/module commands recommended by the local administrators.

No RCP input keyword depends on sham. The same physics input can be executed on a supported workstation or a different cluster after CP2K is built correctly there.
