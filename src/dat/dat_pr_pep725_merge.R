#'******************************************************************************
#' Description: Merge species, station, elevation information of PEP725 dataset 
#' into one gaint table.
#'******************************************************************************

library(data.table)

# in:
# PEP725 stations
pep_stations_file <- file.path(hpc_local_dir, "Data/ARD/pep725_station.csv")
# PEP725 station elevation
pep_elev_file <- file.path(
    hpc_local_dir, 
    "Data/ARD/pep725_station_elevation.csv"
)

# out: pep725_eobs.csv
pep_obs_file <- file.path(gdir, "Data/ARD", "pep725_eobs.csv")

# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Species name
spec_names <- c("aesculus", "alnus", "betula", "fagus", "fraxinus", "quercus")

# Station
pep_stations <- fread(pep_stations_file)
# Elevation
pep_elev <- fread(pep_elev_file)


# Do each species
final_dt <- data.table()
for (i in spec_names) {
    # i <- spec_names[1]
    pep_dt <- fread(file.path(
        hpc_local_dir, 
        "Data/ARD/pep725_pheno_temperature", 
        paste0(i, ".csv"))
    )
    pep_dt$Species <- i
    pep_dt <- merge(pep_dt, pep_elev, by = "PEP_ID")
    pep_dt <- merge(pep_dt, pep_stations, by = "PEP_ID")
    pep_dt[, Tcor := -0.64 * (ALT - ELEV) / 100 + Tmean]

    sel_dt <- NULL
    for (j in unique(pep_dt$PEP_ID)) {
        # j <- 1
        if (length(unique(pep_dt[PEP_ID == j, PhenoYear])) < 40) {
            next
        }
        sel_dt <- rbind(sel_dt, pep_dt[PEP_ID == j, ])
    }
    final_dt <- rbind(final_dt, sel_dt)
}

# Export
fwrite(final_dt, pep_obs_file)

