#!/bin/bash
#SBATCH -J hfhb_fit
#SBATCH -o pipe/hpc_log/hfhb_fit.out
#SBATCH -e pipe/hpc_log/hfhb_fit.err
#SBATCH -N 1
#SBATCH -c 36
#SBATCH -t 0-36:00
#SBATCH -p shared
#SBATCH --mem=16000

srun -c 36 singularity exec ${HOME}/j_r_env_container.sif Rscript src/mod_compare_hf_hb.R

