#!/bin/bash
#SBATCH -J hfhb_cv
#SBATCH -o pipe/hpc_log/hfhb_cv.out
#SBATCH -e pipe/hpc_log/hfhb_cv.err
#SBATCH -N 1
#SBATCH -c 36
#SBATCH -t 0-48:00
#SBATCH -p shared
#SBATCH --mem=16000

srun -c 36 singularity exec ${HOME}/j_r_env_container.sif Rscript src/mod_compare_hf_hb_cv.R

