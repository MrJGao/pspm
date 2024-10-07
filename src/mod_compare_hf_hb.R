#'******************************************************************************
#' Description: Model comparison for ground observed phenology and temperature 
#' at Havard forest and Hubbard Brook Experimental Forest.
#'******************************************************************************
source("src/mod_compare_base.R")

library(parallel)


# ~ Pooled #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HF
hf_pheno_dt <- fread("data/hf_grd_temperature_pheno.csv")
hf_pheno_dt[, Date := as_date(Date)]

# HB
hb_pheno_dt <- fread("data/hb_grd_temperature_pheno.csv")
hb_pheno_dt[, Date := as_date(Date)]

hb_hf_dt <- rbind(
    hf_pheno_dt[, .(siteID, Date,
        Tmean = Temp, Tmax = NA, Tmin = NA, SOS,
        PhenoYear
    )],
    hb_pheno_dt[, .(siteID, Date, Tmean, Tmax, Tmin, SOS, PhenoYear)]
)

print("Fit pooled models....")

# Format data
hb_hf_li <- FormatDataForPhenoModel(hb_hf_dt)
# out: save for later use
saveRDS(hb_hf_li, "pipe/hb_hf_li.Rds")

mod_hb_hf <- FitCompareModels(hb_hf_li)
# out: model fit for ground_pheno + ground_temperature at HF & HB
saveRDS(mod_hb_hf, "pipe/goodness-of-fit/mod_hb_hf_pooled.Rds")




# ~ Site-specific #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print("Fit site-specific models....")

# Make cluster and fit site-specific models
# cl <- makeCluster(detectCores())
# calls <- clusterCall(cl, function() {
#     suppressWarnings({
#         source("src/mod_compare_base.R")
#     })
# })
# clusterExport(cl, c("hb_hf_dt"))
# mod_hb_hf <- clusterApply(cl = cl, unique(hb_hf_dt$siteID), function(x) {
#     cur_site <- hb_hf_dt[siteID == x, ]
#     if (length(unique(cur_site$PhenoYear)) < 15) {
#         return(NULL)
#     }

#     cur_li <- FormatDataForPhenoModel(cur_site)
#     mod <- FitCompareModels(cur_li)

#     return(mod)
# })
# names(mod_hb_hf) <- unique(hb_hf_dt$siteID)


mod_hb_hf <- lapply(unique(hb_hf_dt$siteID), function(x) {
    cur_site <- hb_hf_dt[siteID == x, ]
    if (length(unique(cur_site$PhenoYear)) < 15) {
        return(NULL)
    }

    cur_li <- FormatDataForPhenoModel(cur_site)
    mod <- FitCompareModels(cur_li)
cat(x, "/n")
    return(mod)
})
names(mod_hb_hf) <- unique(hb_hf_dt$siteID)

# out: model fit for MCD12Q2 + Daymet at flux sites
saveRDS(mod_hb_hf, "pipe/goodness-of-fit/mod_hb_hf_site.Rds")







