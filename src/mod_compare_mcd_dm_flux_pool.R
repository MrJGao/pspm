#'******************************************************************************
#' Description: Model comparion for MCD12Q2 + Daymet at flux sites.
#'******************************************************************************

source("src/mod_compare_base.R")


# ~ Pooled ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read the saved Rds
flux_li <- readRDS("pipe/flux_li.Rds")

# Fit pooled models
mod_mcd_dm_flux_pooled <- FitCompareModels(flux_li)

# out: model fit for MCD12Q2 + Daymet at flux sites
saveRDS(mod_mcd_dm_flux_pooled, 
    file.path("pipe/goodness-of-fit", "mod_mcd_dm_flux_pooled",
        "mod_mcd_dm_flux_pooled.Rds"
    )
)

