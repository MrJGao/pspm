#!/bin/bash
#BSUB -J mod_cv
#SBATCH -o pipe/hpc_log/mcd_dm_cv_site.out
#SBATCH -e pipe/hpc_log/mcd_dm_cv_site.err
#SBATCH -N 1
#SBATCH -c 46
#SBATCH -t 0-24:00
#SBATCH -p shared
#SBATCH --mem=16000

srun -c 46 singularity exec ${HOME}/j_r_env_container.sif Rscript src/mod_compare_mcd_dm_flux_cv_site.R ${SLURM_ARRAY_TASK_ID}

