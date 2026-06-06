  GNU nano 5.6.1                                                                                                  slurm/run.slurm                                                                                                             
#!/bin/bash

#SBATCH --job-name=my_test_job
#SBATCH --output=slurm/outputs/slurm-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=aoraki




echo "Sparsh's job started on $(hostname) at $(date)"
echo "Working directory: $(pwd)"

julia scripts/run.jl

echo "Job finished at $(date)"
