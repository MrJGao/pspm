#'******************************************************************************
#' Download spring phneology data and Daymet temperature data for flux sites.
#'******************************************************************************
source("src/base.R")

source("src/hlp/AppEEARS_api.R")


# in: 
# Flux site locations
flux_site_file <- "data/raw/sites_fluxnet2015_Tier1.csv"

# out: data/Fluxsites/
flux_dl_dir <- "data/Fluxsites"



# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read in flux site locations
fluxsites <- fread(flux_site_file)
# Format the location data frame for uploading
pts_df <- fluxsites[, .(
    id = siteID, longitude = lon, latitude = lat,
    category = IGBP
)]


# Download MCD12Q2 and Daymet data using AppEEARS
account <- ParsePrivate(attr_name = "appeears_account")
token <- Login(usr = account$username, pwd = account$password)

# For now, download MCD12Q2's MidGreenup, Daymet, and MCD12Q1 for the IGBP
layers <- data.frame(
    product = c("MCD12Q2.006", rep("DAYMET.004", 2), "MCD12Q1.006"),
    layer = c("MidGreenup", "tmax", "tmin", "LC_Type1")
)

# Submit the job
# Submit a one big job will exceed the max limit of AppEEARS, so split the job.
SubmitTask(token, "Pheno_Daymet_for_Fluxsites",
    task_type = "point", start_date = "2000-01-01", end_date = "2010-12-31",
    layers = layers, point_df = pts_df
)
SubmitTask(token, "Pheno_Daymet_for_Fluxsites_2",
    task_type = "point", start_date = "2011-01-01", end_date = "2020-12-31",
    layers = layers, point_df = pts_df
)

CheckTaskStatus(token, limit = 2, brief = TRUE)


# Download data
DownloadTask(token, task_id = "29ea7a0a-d32d-4bdf-81b7-6ac644c8afd8", 
    out_dir = flux_dl_dir)
DownloadTask(token, task_id = "647231bf-2002-4f19-b6ef-a0ea2241404e", 
    out_dir = flux_dl_dir)
