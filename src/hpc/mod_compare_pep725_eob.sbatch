#!/bin/bash
#BSUB -J pep_al
#SBATCH -o pipe/hpc_log/pep_al.out
#SBATCH -e pipe/hpc_log/pep_al.err
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH -p shared
#SBATCH --mem=4000

srun -c 1 singularity exec ${HOME}/j_r_env_container.sif Rscript src/mod_compare_pep725_eobs.R ${SLURM_ARRAY_TASK_ID}

