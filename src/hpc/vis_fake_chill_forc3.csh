#!/bin/tcsh
if (! -d pipe/log/) then
    mkdir pipe/log/
endif

#BSUB -n 1
#BSUB -W 2:00
#BSUB -R select[avx2]
#BSUB -R "rusage[mem=32]"
#BSUB -J viscf[1-100]
#BSUB -o pipe/log/viscf_out
#BSUB -e pipe/log/viscf_err

conda activate /usr/local/usrapps/jmgray2/jgao/my_env

Rscript src/vis/vis_fake_chill_forc3.R $LSB_JOBINDEX

conda deactivate
