#!/bin/tcsh
if (! -d pipe/hpc_job_log/) then
    mkdir pipe/hpc_job_log/
endif

#BSUB -n 100
#BSUB -W 8:00
#BSUB -J cl_pep
#BSUB -R select[avx2]
#BSUB -oo pipe/hpc_job_log/cl_pep_out
#BSUB -eo pipe/hpc_job_log/cl_pep_err

module load openmpi-gcc/openmpi1.8.4-gcc4.8.2
conda activate my_env

mpirun -n 1 Rscript src/dat/dat_cl_hpc_pep725_temperature.R

conda deactivate
