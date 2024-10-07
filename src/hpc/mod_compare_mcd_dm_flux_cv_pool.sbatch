#!/bin/bash
#BSUB -J mod_cv_pool
#SBATCH -o pipe/hpc_log/mcd_dm_cv_pool.out
#SBATCH -e pipe/hpc_log/mcd_dm_cv_pool.err
#SBATCH -N 1
#SBATCH -c 47
#SBATCH -t 0-7:00
#SBATCH -p shared
#SBATCH --mem=16000

srun -c 47 singularity exec ${HOME}/j_r_env_container.sif Rscript src/mod_compare_mcd_dm_flux_cv_pool.R ${SLURM_ARRAY_TASK_ID}

