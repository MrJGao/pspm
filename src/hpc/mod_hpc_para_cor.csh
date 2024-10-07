#!/bin/tcsh
if (! -d pipe/hpc_job_log/) then
    mkdir pipe/hpc_job_log/
endif

#BSUB -n 102
#BSUB -W 72:00
#BSUB -J par6
#BSUB -oo pipe/hpc_job_log/mod_para_cor_out
#BSUB -eo pipe/hpc_job_log/mod_para_cor_err

module load openmpi-gcc/openmpi1.8.4-gcc4.8.2
conda activate my_env

mpirun -n 1 Rscript src/mod_hpc_para_cor.R

conda deactivate
