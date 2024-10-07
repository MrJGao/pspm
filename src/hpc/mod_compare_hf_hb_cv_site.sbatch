#!/bin/bash
#BSUB -J hfhb_cv_site
#SBATCH -o pipe/hpc_log/hfhb_cv_site.out
#SBATCH -e pipe/hpc_log/hfhb_cv_site.err
#SBATCH -N 1
#SBATCH -c 47
#SBATCH -t 0-36:00
#SBATCH -p shared
#SBATCH --mem=6000

srun -c 47 singularity exec ${HOME}/j_r_env_container.sif Rscript src/mod_compare_hf_hb_cv_site.R ${SLURM_ARRAY_TASK_ID}

