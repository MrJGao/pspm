#'******************************************************************************
#' Description: Survival analysis on MCD12Q2 phenology w/ Daymet and ground 
#' measured temperature.
#'******************************************************************************
rm(list=ls())

# bsub -n 1 -W 24:00 -o Pipeline/hpc_job_log/hfhb_chill_out -e Pipeline/hpc_job_log/hfhb_chill_err -J "chill_mcd2" "Rscript Code/mod_surv_mcd.R"

source("src/base.R")
source("src/hlp/hlp_chilling_detection.R")


# ~ MCD12Q2 + Daymet ####
# ~ ----------------------------------------------------------------------------
mcd_dm_dt <- fread(file.path(hpc_dir, "Data/ARD", 
    "flux_mcd12q2_daymet_dt.csv")
)
mcd_dm_dt[, Date := as.Date(Date)]


# ~ Pooled analysis ----------------------------------------
pool_chill_res <- DetectChillFix(mcd_dm_dt)

# ~ Site-specific analysis ----------------------------------------
sp_chill_res <- lapply(unique(mcd_dm_dt$siteID), function (site_id) {
    site_dt <- mcd_dm_dt[siteID == site_id]
    site_chill_res <- DetectChillFix(site_dt)

    return(site_chill_res)
})


# Combine result
res <- list(
    pool_chill_res = pool_chill_res,
    sp_chill_res = sp_chill_res
)


saveRDS(res, file.path(hpc_dir, "Pipeline", "mcd_dm_chill_res.Rds"))



# ~ MCD12Q2 + flux temperature ####
# ~ ----------------------------------------------------------------------------
mcd_grd_li <- readRDS(file.path(hpc_dir, "Pipeline", "flux_temp_li.Rds"))
mcd_grd_dt <- PhenoModelDataToDT(mcd_grd_li)
mcd_grd_dt[, Date := as.Date(Date)]

# ~ Pooled analysis  ----------------------------------------
pool_chill_res <- DetectChillFix(mcd_grd_dt)

# ~ Site-specific analysis ----------------------------------------
sp_chill_res <- lapply(unique(mcd_grd_dt$siteID), function (site_id) {
    site_dt <- mcd_grd_dt[siteID == site_id]
    site_chill_res <- DetectChillFix(site_dt)

    return(site_chill_res)
})


# Combine result
res <- list(
    pool_chill_res = pool_chill_res,
    sp_chill_res = sp_chill_res
)

saveRDS(res, file.path(hpc_dir, "Pipeline", "mcd_grd_chill_res.Rds"))


