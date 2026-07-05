#!/bin/bash

#SBATCH --job-name=echo-sweep
#SBATCH --output=data/slurm/echo-sweep-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=aoraki
#SBATCH --mail-user=schandra2021nz@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:?Submit this job from the project root}"
mkdir -p data/slurm data/echo_sweeps

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export GKSwstype=100

echo "Echo sweep started on $(hostname) at $(date)"
echo "Working directory: $(pwd)"

julia --project=. scripts/run_echo_parameter_sweep.jl

echo "Echo sweep finished at $(date)"
