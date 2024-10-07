#!/bin/tcsh
if (! -d pipe/hpc_job_log/) then
    mkdir pipe/hpc_job_log/
endif

conda activate my_env


# Goodness of fit
echo "Submit job for goodness-of-fit"
bsub -n 6 -W 3:00 -oo pipe/hpc_job_log/spec_sim_fit_out -eo pipe/hpc_job_log/spec_sim_fit_err "Rscript src/mod_hpc_pspm_spec_sim_fit.R"


# Cross-validation
set num_spec = 5
set i = 1

while ($i <= $num_spec) 
    echo "Submit job index = $i"

    bsub -n 10 -W 3:00 -R "span[hosts=1]" -oo pipe/hpc_job_log/spec_sim_cv_out -eo pipe/hpc_job_log/spec_sim_cv_err "Rscript src/mod_hpc_pspm_spec_sim_cv.R $i"

    @ i++
end



# Test models on mixed data
echo "Submit job for mixed data"
set num_mix_run = 100
set run_idx = 1

while ($run_idx <= $num_mix_run)
    echo "Submit job for mixed data index = $run_idx"

    bsub -n 16 -W 20:00 -R "span[hosts=1]" -oo pipe/hpc_job_log/spec_mix_out -eo pipe/hpc_job_log/spec_mix_err "Rscript src/mod_hpc_pspm_spec_sim_mix.R $run_idx"

    @ run_idx++
end



conda deactivate