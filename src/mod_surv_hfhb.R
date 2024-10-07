#'******************************************************************************
#' Description: Survival analysis on HF & HB dataset.
#'******************************************************************************

# bsub -n 1 -W 24:00 -o pipe/hpc_job_log/hfhb_chill_out -e pipe/hpc_job_log/hfhb_chill_err -J "chill_hfhb" "Rscript Code/mod_surv_hfhb.R"

rm(list=ls())

source("src/base.R")
source("src/hlp/hlp_chilling_detection.R")


hfhb_dt <- fread(file.path(hpc_dir, "Data/ARD", "hb_hf_pheno_temp.csv"))
hfhb_dt[, Date := as.Date(Date)]


# ~ Pooled analysis ####
# ~ ----------------------------------------------------------------------------
pool_chill_res <- DetectChillFix(hfhb_dt)


# ~ Site-specific analysis ####
# ~ ----------------------------------------------------------------------------
sp_chill_res <- lapply(unique(hfhb_dt$siteID), function (site_id) {
    site_dt <- hfhb_dt[siteID == site_id]
    site_chill_res <- DetectChillFix(site_dt)

    return(site_chill_res)
})


# Combine result
res <- list(
    pool_chill_res = pool_chill_res,
    sp_chill_res = sp_chill_res
)

saveRDS(res, file.path(hpc_dir, "Pipeline", "hfhb_chill_res.Rds"))



