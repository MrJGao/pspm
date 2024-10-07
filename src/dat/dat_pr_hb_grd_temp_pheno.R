#'******************************************************************************
#' Clean ground measured air temperature for Hubbard Brook Experimental Forest
#' downloaded from their website, add ground observed phenology preprocessed by
#' the BLSP project.
#'******************************************************************************

# The data zipfile was downloaded from 
# https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-hbr&identifier=59

source("src/base.R")

library(data.table)
library(xml2)
library(sp)
library(raster)
library(lubridate)


# in: 
# Air temperature data at HB directory
hb_temp_data_dir <- "data/raw/HB/air_temperature"
# HB temperature site metadata
hb_temp_data_file <- file.path(
    hb_temp_data_dir,
    "knb-lter-hbr.59.10", 
    "knb-lter-hbr.59.10.xml"
)
# HB temperature data
hb_temp_file <- file.path(
    hb_temp_data_dir,
    "knb-lter-hbr.59.10/HBEF_air_temp_daily_1957-2021.csv"
)
# HF & HB pheno_daymet info
hb_pheno_sites_file <- "data/hf_hb_pheno_daymet.csv"
# Ground phenology at HB sites, from the BLSP project
hb_pheno_file <- "Q:/My Drive/Research/Landsat_LSP/Pipeline/phn_phenos.Rds"


# out: hb_grd_temperature_pheno.csv
grd_temp_pheno_file <- "data/hb_grd_temperature_pheno.csv"


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Check sites to use ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract the site locations from metadata to check if these sites are consistent 
# with phenological observatoinal sites.
hb_meta <- read_xml(hb_temp_data_file)
coverage_node <- xml_find_all(hb_meta, ".//geographicCoverage")
hb_temp_sites <- NULL
for (i in 2:length(coverage_node)) {
    node <- xml_children(coverage_node[i])
    siteID <- xml_text(node[1])
    lon <- xml_text(xml_children(node[2]))[1]
    lat <- xml_text(xml_children(node[2]))[3]
    hb_temp_sites <- rbind(hb_temp_sites, data.frame(
        siteID = siteID, lon = lon, lat = lat))
}
hb_temp_sites <- data.table(hb_temp_sites)
hb_temp_sites[, ":="(lon = as.numeric(as.character(lon)), 
    lat = as.numeric(as.character(lat)))]

hb_temp_sites_sp <- SpatialPointsDataFrame(
    coords = hb_temp_sites[, c("lon", "lat")], 
    proj4string = CRS("+init=epsg:4326"), 
    data = hb_temp_sites
)


# Phenology sites
hb_pheno_sites <- fread(hb_pheno_sites_file)
hb_pheno_sites <- unique(
    hb_pheno_sites[siteID != "US-Ha1", .(siteID, Lon, Lat)]
)
hb_pheno_sites_sp <- SpatialPointsDataFrame(
    coords = hb_pheno_sites[, c("Lon", "Lat")], 
    proj4string = CRS("+init=epsg:4326"), 
    data = hb_pheno_sites
)


# mapview::mapview(list(hb_temp_sites_sp, hb_pheno_sites_sp), 
#     col.regions = c("red", "blue"))

# Ok, the sites are not exactly the same, so I'll use a buffer to extract the 
# sites that are close to the phenology sites
hb_pheno_sites_sp_buf <- buffer(hb_pheno_sites_sp, 500, dissolve = FALSE)

# mapview::mapview(
#     list(hb_pheno_sites_sp_buf, 
#         hb_temp_sites_sp, hb_pheno_sites_sp
#     ), 
#     col.regions = c("blue", "red", "blue")
# )

inter_sp <- raster::intersect(hb_temp_sites_sp, hb_pheno_sites_sp_buf)

# So, based on the intersct and the map, I choose pairs:
# [STA_1, 1B], [STA_6, 5T], [STA14, 7T], [STA_HQ, HQ]


# ~ Extract site temperature data ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hb_temp_dt <- fread(hb_temp_file)
hb_temp_dt <- hb_temp_dt[STA %in% c("HQ", "STA1", "STA6", "STA14")]

# Ground phenology from Bayesian LSP project
hb_pheno <- readRDS(hb_pheno_file)
# Only need midgup
hb_pheno <- hb_pheno[, .(siteID = site, SOS = yday(midgup), 
    PhenoYear = year(midgup))]
hb_pheno <- na.omit(hb_pheno)
# Correct siteID name
hb_pheno$siteID <- gsub("X", "", hb_pheno$siteID)

# Only need 4 sites
hb_pheno <- hb_pheno[siteID %in% c("1B", "5T", "7T", "HQ")]

# Correct siteID
hb_temp_dt$siteID = lapply(hb_temp_dt$STA, function(x) {
    switch(x,
        "HQ" = "HQ",
        "STA1" = "1B",
        "STA6" = "5T",
        "STA14" = "7T"
    )
})
hb_temp_dt[, date := as_date(date)]

# Format temperature and phenology
grd_pheno_dt <- data.table()
for (i in 1:nrow(hb_pheno)) {
    cur_sy <- hb_pheno[i, ]
    # NOTE: The temperature calendar is based on a standard calendar year.
    # All years, including leap years, have 1 - 365 days.

    # R uses a 0 based index for dates
    end_date <- as_date(
        255 - 1, 
        origin = paste0(cur_sy[["PhenoYear"]], "-01-01")
    )
    start_date <- end_date - 364
    # if (leap_year(as_date(paste0(as.numeric(cur_sy[["PhenoYear"]]) - 1, 
    #     "-01-01")))) {
    #     start_date <- start_date - 1
    # }
    cur_temp <- hb_temp_dt[siteID == cur_sy$siteID & 
        between(date, start_date, end_date), ]
    cur_temp$SOS <- cur_sy$SOS
    cur_temp$PhenoYear <- cur_sy$PhenoYear

    grd_pheno_dt <- rbind(grd_pheno_dt, cur_temp)
}
grd_pheno_dt <- grd_pheno_dt[, .(siteID, SOS, PhenoYear, Date = date, 
    Tmean = AVE, Tmax = MAX, Tmin = MIN, STA, Flag)]


# Export HB ground temperature and pheno
fwrite(grd_pheno_dt, grd_temp_pheno_file)
