#!/bin/bash
#BSUB -J pep_cv
#SBATCH -o pipe/hpc_log/pep_cv.out
#SBATCH -e pipe/hpc_log/pep_cv.err
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 0-72:00
#SBATCH -p shared
#SBATCH --mem=4000

srun -c 1 singularity exec ${HOME}/j_r_env_container.sif Rscript src/mod_compare_pep725_eobs_cv.R ${SLURM_ARRAY_TASK_ID}

