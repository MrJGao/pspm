#'******************************************************************************
#' Download Daymet temperature data for Harvard Forest and Hubbard Brook
#' Experimental Forest
#'******************************************************************************
library(data.table)
library(sp)
library(raster)
library(xml2)
source("src/base.R")
source("src/hlp/AppEEARS_api.R")


# in:
# Flux site metadata
site_meta_file <- "data/raw/sites_fluxnet2015_Tier1.csv"
# HB phenological observations
hb_pheno_data_dir <- "data/raw/HB_pheno_obs"

# out: HF_HB_daymet/raw/
data_dl_dir <- "data/raw/HF_HB_daymet"


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Point locations ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Flux site meta ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Flux site meta
site_meta <- fread(site_meta_file)
hf_site <- site_meta[siteID == "US-Ha1"]

hf_site_sp <- hf_site
coordinates(hf_site_sp) <- ~ lon + lat
crs(hf_site_sp) <- "+init=epsg:4326"
# mapview::mapview(hf_site_sp)

# ~ Hubbard brook experimental forest sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parse site coordinates from the xml metadata file
hb_meta <- xml2::read_xml(file.path(
    hb_pheno_data_dir,
    "knb-lter-hbr.51.12/knb-lter-hbr.51.12.xml"
))
coverage_node <- xml_find_all(hb_meta, ".//geographicCoverage")
hb_sites <- NULL
for (i in 2:length(coverage_node)) {
    node <- xml_children(coverage_node[i])
    siteID <- strsplit(xml_text(node[1]), " ")[[1]][3]
    lon <- xml_text(xml_children(node[2]))[1]
    lat <- xml_text(xml_children(node[2]))[3]
    hb_sites <- rbind(
        hb_sites, 
        data.frame(siteID = siteID, lon = lon, lat = lat)
    )
}
hb_sites <- data.table(hb_sites)
hb_sites[, ":="(lon = as.numeric(as.character(lon)), 
    lat = as.numeric(as.character(lat)))]
hb_sites_sp <- SpatialPoints(coords = hb_sites[, c("lon", "lat")])
crs(hb_sites_sp) <- "+init=epsg:4326"
# mapview::mapview(hb_sites_sp)


# ~ Make sample pts dataframe ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samp_pts_df <- rbind(hf_site[, .(siteID, lon, lat)], hb_sites)
colnames(samp_pts_df) <- c("id", "longitude", "latitude")

layers <- data.frame(
    product = rep("DAYMET.004", 2),
    layer = c("tmax", "tmin")
)


# Login to AppEEARS
account <- ParsePrivate(attr_name = "appeears_account")
token <- Login(usr = account$username, pwd = account$password)

# Submit the job
# Submit a one big job will exceed the max limit of AppEEARS, so split the job.
SubmitTask(token, "Daymet_HF_HB_1",
    task_type = "point", start_date = "1988-01-01", end_date = "2000-12-31",
    layers = layers, point_df = samp_pts_df
)
SubmitTask(token, "Daymet_HF_HB_2",
    task_type = "point", start_date = "2001-01-01", end_date = "2020-12-31",
    layers = layers, point_df = samp_pts_df
)


# ~ Check and download #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(jobs <- CheckTaskStatus(token, limit = 3, brief = TRUE))

# Download data
DownloadTask(token, 
    task_id = jobs$task_id[grep("Daymet_HF_HB_1", jobs$task_name)], 
    out_dir = data_dl_dir
)
DownloadTask(token, 
    task_id = jobs$task_id[grep("Daymet_HF_HB_2", jobs$task_name)], 
    out_dir = data_dl_dir
)
