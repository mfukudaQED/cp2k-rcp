#!/bin/bash
#SBATCH -J cp2k_rcp
#SBATCH -p sham
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --cpus-per-task=1
#SBATCH -t 00:10:00

# Template for the ISSP sham cluster. Run as:
# sbatch --export=ALL,CP2K_ENV=/absolute/path/to/cp2k/install/cp2k_env \
#   /path/to/run_cp2k_sham.slurm input.inp output.out
set -eo pipefail

if [ -z "$CP2K_ENV" ] || [ ! -f "$CP2K_ENV" ]; then
  echo "ERROR: set CP2K_ENV to an existing, absolute cp2k_env path." >&2
  exit 2
fi

input="$1"
output="$2"

module purge
source "$CP2K_ENV"

ulimit -l unlimited
ulimit -s unlimited

export I_MPI_COLL_EXTERNAL=no
export I_MPI_FABRICS=shm:ofi
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

mpiexec -n "$SLURM_NTASKS" cp2k.psmp -i "$input" -o "$output"
