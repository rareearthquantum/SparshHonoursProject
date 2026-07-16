#!/bin/bash

#SBATCH --job-name=echo-2d-propagation
#SBATCH --output=data/slurm/echo-2d-propagation-%j.out
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
mkdir -p data/slurm_data/echo_2d_propagation

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export GKSwstype=100

echo "Echo sweep started on $(hostname) at $(date)"
echo "Working directory: $(pwd)"

julia --project=. scripts/run_2d_echo_propagation.jl

echo "Echo 2D propagation finished at $(date)"
