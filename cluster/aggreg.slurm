#!/bin/bash

#SBATCH --job-name=aggreg            # job name
#SBATCH --ntasks=1                        # ask for 1 task
#SBATCH --partition=bigmem
#SBATCH --mem=4G
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --time=02:30:00                   # time limit in hrs:min:sec
#SBATCH --output=test_%j.out       # output and error log


SCRIPT="/home/duchenne/pheno/extract_network_indices.r"
#singularity exec /home/duchenne/ubuntu.sif julia $SCRIPT --args $SLURM_ARRAY_TASK_ID
singularity exec /home/duchenne/ubuntu2.sif R --vanilla -f $SCRIPT --args $SLURM_ARRAY_TASK_ID
#singularity exec /home/duchenne/container_julia.sif julia $SCRIPT $SLURM_ARRAY_TASK_ID