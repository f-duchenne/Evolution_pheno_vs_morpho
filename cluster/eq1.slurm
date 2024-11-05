#!/bin/bash

#SBATCH --job-name=sim_eq1            # job name
#SBATCH --ntasks=1                        # ask for 1 task
#SBATCH --mem=5G 
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00                   # time limit in hrs:min:sec
#SBATCH --output=test_%j.out       # output and error log


SCRIPT="/home/duchenne/pheno/extend_one_trait.jl"
#singularity exec /home/duchenne/ubuntu.sif julia $SCRIPT --args $SLURM_ARRAY_TASK_ID
module load julia-1.6.1-gcc-10.2.0-m3wac3i
sleep $[ $RANDOM % 50 ]
srun julia $SCRIPT $SLURM_ARRAY_TASK_ID
#singularity exec /home/duchenne/container_julia.sif julia $SCRIPT $SLURM_ARRAY_TASK_ID