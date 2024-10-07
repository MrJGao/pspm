#'******************************************************************************
#' Description: Clean and extract the E-OBS temperature data for observational 
#' stations of PEP725 phenology dataset. This script takes a long time, 
#' so I make it run on HPC.
#'******************************************************************************

# Submit the job by the below script
# bsub < Code/dat_cl_hpc_eobs.csh
rm(list=ls())

library(data.table)
library(ncdf4)
library(sp)
library(raster)
library(lubridate)

library(Rmpi)
library(parallel)
library(snow)

source("src/base.R")

wd <- hpc_dir
setwd(wd)

# in:
# Temperature data file
temp_nc <- "Data/EOBS/tg_ens_mean_0.1deg_reg_v24.0e.nc"
# Elevation data file. 
# The v24.0e elevation file on the website is wrong. So, I use v23.1e instead.
# elev_nc <- "Data/EOBS/elev_ens_0.1deg_reg_v24.0e.nc"
elev_nc <- "Data/EOBS/elev_ens_0.1deg_reg_v23.1e.nc"
# PEP725 observations
pep_obs_file <- "Data/ARD/pep725_data.csv"
# PEP725 stations
pep_stations_file <- "Data/ARD/pep725_station.csv"

# out: pep725_eobs_temperature/
# E-OBS temperature at PEP725 stations
out_dir <- "Data/ARD/pep725_eobs_temperature"
# out: pep725_station_elevation.csv
# PEP725 stations elevation
pep_stations_elev_file <- "Data/ARD/pep725_station_elevation.csv"


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PEP725 observations
pep_obs <- fread(pep_obs_file)

# Convert PEP725 stations to spatial vector
pep_stations <- fread(pep_stations_file)
coordinates(pep_stations) <- ~ LON + LAT
proj4string(pep_stations) <- "+init=epsg:4326"


# ~ Extract daily temperature
# To extract temperature records for the observations, I need to find the 
# corresponding observation-year pair to limit the result data as small 
# as possible. Also, since the temperature data cube is chunked by day, 
# to improve the extraction speed, I do the extraction year by year.

temp_nc_in <- nc_open(temp_nc)
temp_lon <- ncvar_get(temp_nc_in, "longitude")
temp_lat <- ncvar_get(temp_nc_in, "latitude")
temp_time <- ncvar_get(temp_nc_in, "time") # time starts from 1950-01-01 00:00
nc_close(temp_nc_in)

temp_year <- data.table(
    temp_time, 
    Date = as_date(temp_time, origin = "1950-01-01")
)
temp_year[, Year := year(Date)]

# Extract temperature
ExtractTemperature <- function(yr) {
    cur_yr <- temp_year[Year == yr, ]
    
    # Subset stations in current year only
    if (yr == 1950) {
        # I need data in 1950 to model phenology of 1951
        cur_pep_obs_id <- pep_obs[YEAR == 1951, PEP_ID]
    } else {
        cur_pep_obs_id <- pep_obs[YEAR == yr, PEP_ID]
    }
    if (length(cur_pep_obs_id) == 0) {
        return(NULL)
    }
    cur_pep_stations <- subset(
        pep_stations, 
        pep_stations$PEP_ID %in% cur_pep_obs_id
    )
    if (nrow(cur_pep_stations@data) == 0) {
        return(NULL)
    }

    temp_nc_in <- nc_open(temp_nc)
    temp_dt <- lapply(cur_yr$temp_time, function(t) {
        # Get curernt day's temperature raster
        tg <- ncvar_get(temp_nc_in, "tg", start = c(1, 1, t + 1), 
            count = c(-1, -1, 1))
        temp_img <- raster(tg,
            xmn = min(temp_lon), xmx = max(temp_lon),
            ymn = min(temp_lat), ymx = max(temp_lat),
            crs = "+init=epsg:4326"
        )

        # Extract values
        val <- extract(temp_img, cur_pep_stations)
        exa_dt <- data.table(
            PEP_ID = cur_pep_stations$PEP_ID,
            Date = cur_yr[temp_time == t, Date],
            Tmean = val
        )
        return(exa_dt)
    })
    nc_close(temp_nc_in)

    temp_dt <- do.call(rbind, temp_dt)

    # Order by ID and convert `time` to Date strings
    setorder(temp_dt, PEP_ID)

    # E-OBS temperature at PEP725 stations
    # As the result data file is too large, it'll affect the speed of later 
    # querying and processing. I split the records by yearly chunks.
    if (!dir.exists(out_dir)) {
        dir.create(out_dir)
    }
    fwrite(temp_dt, file = file.path(out_dir, paste0(yr, ".csv")))
}



# make cluster
cl <- makeCluster((mpi.universe.size() - 1), type = "MPI")
# cl <- makeCluster((23 - 1), type = "SOCK")

# export variables to all nodes
clusterExport(cl, list = c(
    "wd", "out_dir", 
    "temp_nc", "temp_lon", "temp_lat", "temp_year",
    "pep_stations", "pep_obs"
))
# load libraries and set up working directory on all nodes
calls <- clusterCall(cl, function() {
    suppressMessages({
        library(data.table)
        library(ncdf4)
        library(sp)
        library(raster)
        library(lubridate)
    })
    setwd(wd)
})


# cluster run
outputs <- clusterApplyLB(
    cl = cl, x = unique(temp_year[, Year]),
    fun = ExtractTemperature
)


# ~ Extract elevation
elev_nc_in <- nc_open(elev_nc)
elev_lon <- ncvar_get(elev_nc_in, "longitude")
elev_lat <- ncvar_get(elev_nc_in, "latitude")
elev <- ncvar_get(elev_nc_in, "elevation")
elev_img <- raster(elev,
    xmn = min(elev_lon), xmx = max(elev_lon), ymn = min(elev_lat), 
    ymx = max(elev_lat),
    crs = "+init=epsg:4326"
)
elev_val <- extract(elev_img, pep_stations)
elev_dt <- data.table(PEP_ID = pep_stations$PEP_ID, ELEV = elev_val)

nc_close(elev_nc_in)

# Export PEP725 stations' elevation
fwrite(elev_dt, file = pep_stations_elev_file)


# shut down cluster
snow::stopCluster(cl)
Rmpi::mpi.quit()
