#!/bin/tcsh
if (! -d Pipeline/hpc_job_log/) then
    mkdir Pipeline/hpc_job_log/
endif

#BSUB -n 73
#BSUB -W 5:00
#BSUB -J cl_eobs
#BSUB -oo Pipeline/hpc_job_log/cl_eobs_out
#BSUB -eo Pipeline/hpc_job_log/cl_eobs_err

module load openmpi-gcc/openmpi1.8.4-gcc4.8.2
conda activate my_env

mpirun -n 1 Rscript src/dat/dat_cl_hpc_eobs.R

conda deactivate
