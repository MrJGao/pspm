#'******************************************************************************
#' Description: Fit classic spring phenology models to various data source and 
#' save the result to a file that will be used later to make figures.
#'******************************************************************************

# A single run takes about 8 hours b/c the UN model requires 50 runs to find the
# best and each run takes around 9 min.

rm(list=ls())

# library(parallel)
# library(snow)

source("src/mod_compare_base.R")

args = commandArgs(trailingOnly = TRUE)
idx <- as.numeric(args[1])

# Set up cluster
# cl <- makeCluster(ncores)

# calls <- clusterCall(cl, function() {
#     suppressWarnings({
#         source("src/mod_compare_base.R")
#         library(parallel)
#         library(snow)
#     })
# })


# Station
# pep_stations <- fread("data/raw/pep725_station.csv")
# # Elevation
# pep_elev <- fread("data/raw/pep725_station_elevation.csv")


# Chunk each site and species
ChunkSiteSpecies <- function() {
    print("Caching...")

    species_files <- list.files(
        "data/raw/pep725_pheno_temperature",
        ".csv$",
        full.names = TRUE
    )

    if (!dir.exists("pipe/pep_tmp")) {
        dir.create("pipe/pep_tmp")
    }
    for (i in species_files) {
        # i <- species_files[1]
        species_name <- strsplit(basename(i), "\\.")[[1]][1]

        # Calibrate elevation
        pep_dt <- fread(i)
        pep_dt <- merge(pep_dt, pep_elev, by = "PEP_ID")
        pep_dt <- merge(pep_dt, pep_stations, by = "PEP_ID")
        pep_dt[, Tcor := -0.64 * (ALT - ELEV) / 100 + Tmean]
        pep_dt <- pep_dt[, .(
            siteID = PEP_ID, Date, Tmax = NA, Tmin = NA, Tmean = Tcor,
            SOS, PhenoYear
        )]

        for (j in unique(pep_dt$siteID)) {
            # j <- 1
            if (length(unique(pep_dt[siteID == j, PhenoYear])) < 40) {
                next
            }
            saveRDS(pep_dt[siteID == j, ], file.path(
                hpc_dir, "pipe/pep_tmp",
                paste0(species_name, "_", j, ".Rds")
            ))
        }
    }
}

# ChunkSiteSpecies()


# species_name <- "aesculus"
# species_name <- "alnus"
# species_name <- "betula"
# species_name <- "fagus"
# species_name <- "fraxinus"
# species_name <- "quercus"


# Fit species-specific models
# caches <- list.files("pipe/pep_tmp",
#     paste0(species_name),
#     full.names = TRUE
# )


caches <- list.files("pipe/pep_tmp", full.names = TRUE)

outdir <- file.path("pipe/goodness-of-fit", "pep")
if (!dir.exists(outdir)) {
    dir.create(outdir)
}

print(paste("Fit models......", basename(caches[idx])))

cur_site <- readRDS(caches[idx])

cur_li <- FormatDataForPhenoModel(cur_site)
system.time({
    mod <- FitCompareModels(cur_li, if_un_single = TRUE)
})

saveRDS(mod, file.path(outdir, basename(caches[idx])))

# pep_eobs <- data.table()
# pep_eobs <- clusterApplyLB(cl = cl, caches, function(x) {
#     cur_site <- readRDS(x)

#     cur_li <- FormatDataForPhenoModel(cur_site)
#     mod <- FitCompareModels(cur_li, if_un_single = TRUE)

#     return(mod)
# })
# names(pep_eobs) <- lapply(caches, function(x) {
#     strsplit(basename(x), "\\.")[[1]][1]
# })

# out: model fit for PEP725 + E-OBS
# saveRDS(pep_eobs, file.path("pipe/goodness-of-fit", 
#     paste0("mod_", species_name, "_pep725_eobs", ".Rds")))


# # shut down cluster
# snow::stopCluster(cl)
# Rmpi::mpi.quit()
