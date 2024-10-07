#!/bin/tcsh
if (! -d pipe/hpc_job_log/) then
    mkdir pipe/hpc_job_log/
endif

conda activate my_env


set run = 4000
set i = 1

while ($i <= $run) 
    echo "Submit job index = $i"

    bsub -n 1 -W 1:20 -oo pipe/hpc_job_log/sim_noise_out -eo pipe/hpc_job_log/sim_noise_err -J "snoise_$i" "Rscript src/mod_hpc_pspm_sim_noise.R $i"

    @ i++
end


conda deactivate