#'******************************************************************************
#' Description: Download temperature data for Hubbard Brook Experimental Forest.
#'******************************************************************************
library(httr)

source("src/base.R")


# The Data downloading page is:
# https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-hbr&identifier=59
# Users can download it manually, but I'll do it by code.

# in: 
# Data website
url <- "https://portal.edirepository.org/nis/archiveDownload"

# out: HB_air_temperature/
out_dir <- "data/raw/HB_air_temperature"
if (!dir.exists(out_dir)) dir.create(out_dir)


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# HEAD request to get the download filename
res <- HEAD(url,
    query = list(
        packageid = "knb-lter-hbr.59.10",
        archive = "Download Zip Archive"
    )
)
filename <- strsplit(headers(res)$`content-disposition`, "=")[[1]][2]

# GET request to download the data
res <- GET(url,
    query = list(
        packageid = "knb-lter-hbr.59.10",
        archive = "Download Zip Archive"
    ),
    write_disk(file.path(out_dir, filename), overwrite = TRUE)
)

