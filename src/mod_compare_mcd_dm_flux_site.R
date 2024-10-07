#'******************************************************************************
#' Description: Model comparion for MCD12Q2 + Daymet at flux sites.
#'******************************************************************************

source("src/mod_compare_base.R")
library(lubridate)
library(data.table)

args = commandArgs(trailingOnly = TRUE)
idx <- as.numeric(args[1])


# Make cluster and fit site-specific models
# cl <- makeCluster(detectCores() - 1)
# calls <- clusterCall(cl, function() {
#     wd <- getwd()
#     setwd(wd)

#     suppressWarnings({
#         source("Code/mod_compare_base.R")
#     })
# })
# clusterExport(cl, c("flux_dt"))

# Read data in
flux_dt <- fread("data/flux_mcd12q2_daymet_dt.csv")
flux_dt[, Date := as_date(Date)]
# merge site meta
site_meta <- fread("data/raw/sites_fluxnet2015_Tier1.csv")
flux_dt <- merge(flux_dt, site_meta[, .(siteID, IGBP)], by = "siteID")
# filter out SOS > 255
flux_dt <- flux_dt[SOS < 255 & IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), ]

sites <- unique(flux_dt$siteID)
x <- sites[idx]

print(paste("Fit MCD12Q2 + Daymet site-specific....", x))

cur_site <- flux_dt[siteID == x, ]
if (length(unique(cur_site$PhenoYear)) < 10) return(NULL)

cur_li <- FormatDataForPhenoModel(cur_site)
mod <- FitCompareModels(cur_li)

saveRDS(mod, 
    file.path("pipe/goodness-of-fit", "mcd_dm_flux_site", paste0(x, ".Rds"))
)

# mod_mcd_dm_flux <- clusterApply(cl = cl, unique(flux_dt$siteID), function(x) {
#     cur_site <- flux_dt[siteID == x, ]
#     if (length(unique(cur_site$PhenoYear)) < 10) return(NULL)

#     cur_li <- FormatDataForPhenoModel(cur_site)
#     mod <- FitCompareModels(cur_li)

#     return(mod)
# })
# names(mod_mcd_dm_flux) <- unique(flux_dt$siteID)

# # out: model fit for MCD12Q2 + Daymet at flux sites
# saveRDS(mod_mcd_dm_flux, file.path(gdir, "Pipeline/mod_mcd_dm_flux_site.RData"))



