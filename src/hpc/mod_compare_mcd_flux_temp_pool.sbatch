#!/bin/bash
#BSUB -J flux_temp_pool
#SBATCH -o pipe/hpc_log/mcd_flux_temp_pool.out
#SBATCH -e pipe/hpc_log/mcd_flux_temp_pool.err
#SBATCH -N 1
#SBATCH -c 45
#SBATCH -t 0-12:00
#SBATCH -p shared
#SBATCH --mem=16000

srun -c 45 singularity exec ${HOME}/j_r_env_container.sif Rscript src/mod_compare_mcd_flux_temp_pool.R
