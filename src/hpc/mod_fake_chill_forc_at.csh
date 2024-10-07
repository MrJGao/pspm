#!/bin/tcsh
if (! -d pipe/log/) then
    mkdir pipe/log/
endif

#BSUB -n 20
#BSUB -W 12:00
#BSUB -R "rusage[mem=64]"
#BSUB -J atsim
#BSUB -oo pipe/log/atsim_out
#BSUB -eo pipe/log/atsim_err

module load openmpi-gcc/openmpi1.8.4-gcc4.8.2

conda activate /usr/local/usrapps/jmgray2/jgao/my_env
mpirun -n 1 Rscript src/mod_fake_chill_forc_at.R

conda deactivate
