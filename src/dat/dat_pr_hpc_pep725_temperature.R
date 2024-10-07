#'******************************************************************************
#' Description: Format the temperature records for each observations for later 
#' modeling.
#'******************************************************************************

# To run this script, do:
# bsub < src/dat_cl_hpc_pep_pheno_temperature_submit.csh

library(data.table)
library(lubridate)

library(Rmpi)
library(parallel)
library(snow)

source("src/base.R")

wd <- hpc_dir
setwd(wd)

# in:
# PEP725 data
pep_dt_file <- "data/pep725_data.csv"
# PEP725 station elevation
pep_elev_file <- "data/pep725_station_elevation.csv"
# E-OBS temperature at PEP725 stations
temp_data_dir <- "data/pep725_eobs_temperature"

# out: pep725_pheno_temperature/
# PEP725 phenology and temperature
out_dir <- "data/pep725_pheno_temperature"
if (!dir.exists(out_dir)) { dir.create(out_dir) }


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Load and preprocess ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
pep_dt <- fread(pep_dt_file)
pep_elev <- fread(pep_elev_file)
# Merge data with stations
pep_dt <- merge(pep_dt, pep_elev)

# Follow Zhang et al., 2022, I select same 6 species from the data, 
# and I choose `BBCH == 11` as the leaf unfolding date.
pep_dt <- pep_dt[species %in% c(
    "Aesculus hippocastanum",
    "Alnus(Alnus glutinosa)",
    "Betula(Betula pendula_(B._verrucosa__B._alba))",
    "Fagus(Fagus sylvatica)",
    "Fraxinus excelsior",
    "Quercus robur_(Q.peduncula)"
), ]
pep_dt <- pep_dt[BBCH == 11, ]

# Make alias of species names
alias <- vapply(pep_dt$species, function(x) {
    switch(x,
        "Aesculus hippocastanum" = "aesculus",
        "Alnus(Alnus glutinosa)" = "alnus",
        "Betula(Betula pendula_(B._verrucosa__B._alba))" = "betula",
        "Fagus(Fagus sylvatica)" = "fagus",
        "Fraxinus excelsior" = "fraxinus",
        "Quercus robur_(Q.peduncula)" = "quercus"
    )
}, character(1))
pep_dt$alias <- alias

# Extract temperature records for each species and each observation-year
species <- unique(pep_dt$alias)


# ~ Cluster ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make cluster
cl <- makeCluster((mpi.universe.size() - 1), type = "MPI")
# cl <- makeCluster((24 - 1), type = "SOCK")

# export variables to all nodes
clusterExport(cl, list = c(
    "wd", "temp_data_dir", "out_dir"
))

# load libraries and set up working directory on all nodes
calls <- clusterCall(cl, function() {
    suppressMessages({
        library(data.table)
        library(lubridate)
    })
    setwd(wd)
})

for (spec in species) {
    # spec <- species[1]

    spec_dt <- pep_dt[alias == spec, ]

    # export variables to all nodes
    clusterExport(cl, list = c("spec_dt"), envir = environment())

    set.seed(1000)
    samples <- sample(1:nrow(spec_dt), 1e4)

    # cluster run
    outputs <- clusterApplyLB(cl = cl, x = samples, fun = function(i) {
        # obs_yr <- spec_dt[PEP_ID == 2000 & YEAR == 2001]
        obs_yr <- spec_dt[i, ]

        # NOTE: The temperature calendar is based on a standard calendar year.
        # For leap years, the temperature dataset has 366 days.
        # I'm gonna drop the first start counting day for leap years to keep it
        # 365 days.
        # R uses a 0 based index for dates
        end_date <- as_date(
            255 - 1, 
            origin = paste0(obs_yr[["YEAR"]], "-01-01")
        )
        start_date <- end_date - 364

        # Read temperature
        cur_temp <- rbind(
            fread(file.path(temp_data_dir, paste0(year(start_date), ".csv"))),
            fread(file.path(temp_data_dir, paste0(year(end_date), ".csv")))
        )
        cur_temp[, Date := as_date(Date)]

        cur <- cur_temp[PEP_ID == obs_yr$PEP_ID &
            between(Date, start_date, end_date), ]
        # Remove NA
        cur <- cur[!is.na(Tmean), ]
        if (nrow(cur) < 365) { # Temperature data missing
            return(NULL)
        }
        cur$SOS <- obs_yr$DAY
        cur$PhenoYear <- obs_yr$YEAR
        return(cur)
    })
    spec_temp_dt <- do.call(rbind, outputs)

    # Export species specific phenology and temperature
    fwrite(spec_temp_dt, file = file.path(out_dir, paste0(spec, ".csv")))
}


# shut down cluster
snow::stopCluster(cl)
Rmpi::mpi.quit()
