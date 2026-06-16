#!/bin/bash

#SBATCH --job-name=my_test_job
#SBATCH --output=data/slurm/slurm-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=aoraki
#SBATCH --schandra2021nz@gmail.com


mkdir -p data/slurm
mkdir -p plots

echo "Sparsh's job started on $(hostname) at $(date)"
echo "Working directory: $(pwd)"

julia scripts/run.jl

echo "Job finished at $(date)"
