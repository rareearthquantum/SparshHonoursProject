#!/bin/bash

#SBATCH --job-name=my_fft_test_job
#SBATCH --output=slurm/outputs/slurm-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=aoraki
#SBATCH --mail-user=schandra2021nz@gmail.com
#SBATCH --mail-type=END,FAIL



echo "Sparsh's job started on $(hostname) at $(date)"
echo "Working directory: $(pwd)"

julia scripts/prototype_fft_run.jl

echo "Job finished at $(date)"
