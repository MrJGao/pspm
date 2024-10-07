#'******************************************************************************
#' Description: Clean and format PEP725 data.
#'******************************************************************************
rm(list=ls())

library(data.table)

source("src/base.R")

# in: 
# pep725 data directory
data_dir <- file.path(hpc_local_dir, "Data/PEP725")


# out: pep725_data.csv
pep_obs_file <- file.path(gdir, "Data/ARD/pep725_data.csv")
# out: pep725_station.csv
pep_station_file <- file.path(gdir, "Data/ARD/pep725_station.csv")
# out: pep725_bbch.csv
pep_bbch_meta_file <- file.path(gdir, "Data/ARD/pep725_bbch.csv")

# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data files
dat_files <- list.files(data_dir, pattern = "PEP725_.*.tar.gz$")

# Unpack the files to a temporary directory, read the data out, and combine 
# and save the data to my ARD folder
tmp_dir <- tempdir()
dat_dt <- NULL
station_dt <- NULL
for (i in seq_along(dat_files)) {
    tryCatch(
        {
            tarfile <- file.path(data_dir, dat_files[i])

            # To avoid special characters in file names, I untar the entire 
            # file to a temp folder, and find the file names later.
            untar(tarfile, exdir = path.expand(file.path(tmp_dir, "tmp")))
            tarlist <- list.files(path.expand(file.path(tmp_dir, "tmp")))
            data_file <- tarlist[
                -grep(".*._(stations)|(README)|(BBCH).*.", tarlist)
            ]
            station_file <- tarlist[grep(".*._stations.*.", tarlist)]

            # Get the BBCH description, this only needs to be done once
            if (i == 1) {
                bbch_file <- tarlist[grep(".*._BBCH.*.", tarlist)]
                bbch_dt <- fread(file.path(tmp_dir, "tmp", bbch_file))
            }
            # read the data
            dat <- fread(file.path(tmp_dir, "tmp", data_file))
            station <- fread(
                file.path(tmp_dir, "tmp", station_file), 
                fill = TRUE
            )

            # discard any columns > 6 (errors in NAME field)
            if (ncol(station) > 6) {
                station <- station[, -c(7:ncol(station)), with = FALSE]
            }

            # manually assign column names to avoid errors with malformed data
            colnames(station) <- c("PEP_ID", "National_ID", "LON", "LAT", "ALT", 
                "NAME")

            # convert to numeric and add fields
            station[, ":="(LON = as.numeric(LON), LAT = as.numeric(LAT))]

            dat$country <- substr(data_file, 8, 9)
            dat$species <- sub(
                "_", 
                " ", 
                substr(data_file, 11, nchar(data_file) - 4)
            )

            dat_dt <- rbind(dat_dt, dat)
            station_dt <- rbind(station_dt, station)
        },
        error = function(err) {
            print(paste0(dat_files[i], "failed, check it!"))
        },
        finally = {
            # remove the temporary files
            unlink(file.path(tmp_dir, "tmp"), recursive = TRUE)
        }
    )
}

# Stations are repeated
station_dt <- unique(station_dt)

# ! Only one file `PEP725_IE_108_012.tar.gz` failed because of a special 
# ! character in the file name, but the file actually doesn't contain any data, 
# ! so that's fine.


# As the E-OBS temperature is from 1950-2021, we wouldn't be able to use PEP 
# observations out of this time period. So, I only need PEP observations within 
# this time period.
dat_dt <- dat_dt[YEAR > 1950, ]
# And follow Zhang et al 2022, I only need BBCH = 11 as the leaf unfolding date
dat_dt <- dat_dt[BBCH == 11, ]

# Also, I subset stations that are corresponding to these observations
station_dt <- station_dt[PEP_ID %in% dat_dt$PEP_ID, ]


# PEP725 data and corresponding metadata
fwrite(dat_dt, file = pep_obs_file)
fwrite(station_dt, file = pep_station_file)
fwrite(bbch_dt, file = pep_bbch_meta_file)
