#!/bin/tcsh
if (! -d pipe/hpc_job_log/) then
    mkdir pipe/hpc_job_log/
endif

#BSUB -n 100
#BSUB -W 24:00
#BSUB -J iter_sim
#BSUB -oo pipe/hpc_job_log/iter_sim_out
#BSUB -eo pipe/hpc_job_log/iter_sim_err

module load openmpi-gcc/openmpi1.8.4-gcc4.8.2
conda activate my_env

mpirun -n 1 Rscript src/mod_hpc_iteration_sim.R

conda deactivate
