#!/bin/bash
#SBATCH -J rmd
#SBATCH -o 111.out
#SBATCH -e 111.err
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-20:00
#SBATCH -p shared
#SBATCH --mem=16000

# Run the program
srun -c 4 singularity exec ${HOME}/j_r_env_container.sif Rscript -e 'rmarkdown::render("src/test_bayes_fit.Rmd", knit_root_dir = "../")'
