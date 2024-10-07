#'******************************************************************************
#' Description: Maps for phenological data sites.
#'******************************************************************************

source("src/base.R")
source("src/vis/vis_base.R")

library(data.table)
library(lubridate)
library(sp)
library(rnaturalearth)
library(RColorBrewer)
library(grid)
library(tmap)
library(maps)
library(rgeos)
library(xml2)
library(grid)


map_colors <- brewer.pal(8, "Set1")
deci_col <- map_colors[5]
ever_col <- map_colors[3]
mf_col <- map_colors[2]


# ~ MCD12Q2 at forest flux sites ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ With Daymet temperature ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read data in
flux_dt <- fread(file.path(gdir, "Data/ARD/flux_mcd12q2_daymet_dt.csv"))
flux_dt[, Date := as_date(Date)]

# Remove the Alaska site
flux_dt <- flux_dt[!siteID %in% c("US-Prr", "US-Me1", "US-Blo")]

# merge site meta
site_meta <- fread(file.path(gdir, "Data/ARD/sites_fluxnet2015_Tier1.csv"))
flux_dt <- merge(flux_dt, site_meta[, .(siteID, lat, lon, IGBP)], by = "siteID")
# filter out SOS > 255
flux_dt <- flux_dt[SOS < 255 & IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), ]

flux_sites <- unique(flux_dt[, .(siteID, lat, lon, IGBP)])
coordinates(flux_sites) <- ~ lon + lat
proj4string(flux_sites) <- "+init=epsg:4326"


# ~ With flux temperature ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# in: flux temperature
flux_temp <- fread(file.path(
    gdir,
    "Data/ARD/flux_measured_daily_temperature.csv"
))
flux_temp[, Date := as_date(Date)]

# in: flux mcd12q2 phenology
flux_dt <- fread(file.path(gdir, "Data/ARD/flux_mcd12q2_daymet_dt.csv"))
flux_dt[, Date := as_date(Date)]
# Remove the Alaska site
flux_dt <- flux_dt[!siteID %in% c("US-Prr", "US-Me1", "US-Blo")]

# merge site meta
site_meta <- fread(file.path(gdir, "Data/ARD/sites_fluxnet2015_Tier1.csv"))
flux_dt <- merge(flux_dt, site_meta[, .(siteID, IGBP)], by = "siteID")
# filter out SOS > 255
flux_dt <- flux_dt[SOS < 255 & IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), ]
# Merge the two
flux_dt <- merge(flux_dt, flux_temp, by = c("siteID", "Date"), all.x = TRUE)

# Filter out site-years with measured temperature data that have missing values
flux_temp_dt <- by(flux_dt, flux_dt[, .(siteID, PhenoYear)], function(sy) {
    if (any(is.na(sy$Temp)) == FALSE) {
        return(sy)
    }
})
flux_temp_dt <- do.call(rbind, flux_temp_dt)

# Filter out non-forest sites
flux_temp_dt <- flux_temp_dt[IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), .(
    siteID, Date,
    Tmax = NA, Tmin = NA, Tmean = Temp, SOS, PhenoYear
)]



north_america <- ne_countries(continent = "North America", scale = "medium")
bb <- bbox(flux_sites)


# # fig: forest flux site maps
flux_sites_map <- tm_shape(north_america, bbox = c(-150, 25, -65, 60)) +
    tm_polygons() +
    tm_shape(flux_sites) +
    tm_dots(size = 0.1, col = "IGBP", 
        palette = palette(c(deci_col, ever_col, mf_col)), 
        title = "", 
        labels = c("Deciduous forest", "Evergreen forest", "Mixed forest")) +
    tm_shape(subset(flux_sites, siteID %in% unique(flux_temp_dt$siteID))) +
    tm_dots(size = 0.1, shape = 4) + 
    tm_compass(position = c("left", "top")) +
    tm_scale_bar(position = c("left", "bottom")) +
    tm_grid(lines = FALSE) +
    tm_xlab("Longitude") +
    tm_ylab("Latitude") +
    tm_layout(legend.position = c(0.01, 0.2)) +
    tm_add_legend(type = "symbol", 
        labels = "With ground temperature", 
        shape = 4
    )

tmap_save(flux_sites_map, 
    filename = file.path("Output/Assets", "flux_sites_map.png"),
    width = 2600, height = 2000, dpi = 500)




# ~ HF & HB ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grd_pheno_dt <- fread(file.path(gdir, "Data/ARD/hf_hb_pheno_daymet.csv"))
hf_hb_sites <- unique(grd_pheno_dt[, .(siteID, Lon, Lat)])
coordinates(hf_hb_sites) <- ~ Lon + Lat
proj4string(hf_hb_sites) <- CRS("+proj=longlat +datum=WGS84")


# Air temperature sites
hb_meta <- read_xml(file.path(
    gdir, "Data/Raw_data/HB/air_temperature",
    "knb-lter-hbr.59.10", "knb-lter-hbr.59.10.xml"
))
coverage_node <- xml_find_all(hb_meta, ".//geographicCoverage")
hb_temp_sites <- NULL
for (i in 2:length(coverage_node)) {
    node <- xml_children(coverage_node[i])
    siteID <- xml_text(node[1])
    lon <- xml_text(xml_children(node[2]))[1]
    lat <- xml_text(xml_children(node[2]))[3]
    hb_temp_sites <- rbind(hb_temp_sites, data.frame(
        siteID = siteID, lon = lon, lat = lat
    ))
}
hb_temp_sites <- data.table(hb_temp_sites)
hb_temp_sites[, ":="(lon = as.numeric(as.character(lon)),
    lat = as.numeric(as.character(lat)))]

hb_temp_sites_sp <- SpatialPointsDataFrame(
    coords = hb_temp_sites[, c("lon", "lat")],
    proj4string = CRS("+init=epsg:4326"), data = hb_temp_sites
)



states <- ne_states(country = "united states of america")


# fig: HF & HB site map
# HB sites
hb_bbox <- bbox2SP(43.91, 43.97, -71.78, -71.6)
tile <- maptiles::get_tiles(terra::ext(hb_bbox),
    zoom = 12,
    crop = TRUE,
    provider = "Stamen.Terrain"
)


hb_map <- tm_shape(tile, projection = "epsg:2032",
        bbox = c(759000, 4867000, 765000, 4873000)) +
    tm_rgb(alpha = 1) +
    tm_shape(subset(hf_hb_sites, siteID != "US-Ha1")) +
    tm_dots(shape = 17, size = 0.1) +
    tm_shape(hb_temp_sites_sp) +
    tm_dots(shape = 16, col = "blue", size = 0.1) +
    tm_compass(position = c("left", "top")) +
    tm_scale_bar(position = c("left", "bottom"), width = 0.2) +
    tm_grid(lines = FALSE, labels.size = 0.5) +
    tm_xlab("x") +
    tm_ylab("y") + 
    tm_add_legend(type = "symbol", col = "black", size = 0.3, shape = 17,
        labels = "Phenological sites") +
    tm_add_legend(type = "symbol", col = "blue", size = 0.3, shape = 16,
        labels = "Temperature stations") +
    tm_layout(
        title = "Hubbard Brook Experimental Forest", 
        title.position = c(0.27, 0.95),
        title.size = 0.8,
        # frame = mf_col,
        # frame.lwd = 3,
        legend.position = c(0.6, 0.05),
        legend.bg.color = "#dcd8d8")

# HF & HB sites
# hf_hb_bbox <- bbox2SP(42.5, 44, -72.5, -71.1)
# hf_hb_map <- tm_shape(states, bbox = hf_hb_bbox) +
#     tm_polygons(col = "white") +
#     # tm_shape(hb_bbox) +
#     # tm_borders(lwd = 3, col = mf_col) +
#     tm_shape(hf_hb_sites) +
#     tm_dots(shape = 17, size = 0.1) +
#     tm_credits(text = "Harvard Forest", position = c(0.2, 0.05)) +
#     tm_credits(text = "Hubbard Brook\n Experimental Forest", 
#         position = c(0.3, 0.75))
    # tm_layout(frame = ever_col, frame.lwd = 3)

# US boundary
us_base_map <- tm_shape(states, bbox = c(-83, 25, -67, 50)) +
    tm_polygons(border.col = "white", labels = 1) +
    tm_shape(hf_hb_sites) +
    tm_dots(shape = 17, size = 0.1) +
    tm_layout(frame = FALSE) +
    # tm_grid(lines = FALSE, labels.size = 0.5) +
    tm_credits(text = "Harvard Forest", position = c(0.25, 0.63), size = 0.6) +
    tm_credits(text = "Hubbard Brook\nExperimental Forest", 
        position = c(0.063, 0.75), size = 0.6)
    # tm_shape(hf_hb_bbox) +
    # tm_borders(lwd = 3, col = ever_col)



# NOTE: This figure was further processed by drawio
png(file.path("Output/Assets", "hf_hb_map.png"), 
    width = 1800, height = 1350, res = 300
)
# print(hf_hb_map, 
#     vp = viewport(0.13, 0.9, just = "top", width = 0.4, height = 0.4)
# )
print(us_base_map, 
    vp = viewport(0.14, 0.95, just = "top", width = 0.25, height = 1)
)
print(hb_map, 
    vp = viewport(0.62, 0.95, just = "top", width = 0.7, height = 1)
)
grid.lines(
    x = c(0.185, 0.388),
    y = c(0.611, 0.863)
)
grid.lines(
    x = c(0.185, 0.388),
    y = c(0.599, 0.11)
)

dev.off()






# ~ PEP725 stations ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Station
pep_stations <- fread(file.path(gdir, "Data/ARD/pep725_station.csv"))

# Use cached files to select stations with more than 40 site years b/c we only 
# used them to conduct the analysis
cached <- caches <- list.files(
    file.path(
        hpc_local_dir,
        "Pipeline/tmp"
    ),
    ".*.Rds$"
)
# Get unique station IDs
station_ids <- sapply(cached, function(x) {
    strsplit(tools::file_path_sans_ext(basename(x)), "_")[[1]][2]
})
station_ids <- unique(station_ids)

pep_stations <- pep_stations[PEP_ID %in% station_ids]
coordinates(pep_stations) <- ~ LON + LAT
proj4string(pep_stations) <- "+init=epsg:4326"

europe <- ne_countries(continent = "Europe")

# PEP725 stations map
pep725_station_map <- tm_shape(europe, xlim = bbox(pep_stations)[c(1, 3)], 
    ylim = bbox(pep_stations)[c(2, 4)]) +
    tm_polygons(border.col = "white") +
    tm_shape(pep_stations) +
    tm_dots(size = 0.05) +
    tm_compass(position = c("right", "top")) +
    tm_scale_bar(position = c("left", "bottom"), width = 0.3) +
    tm_grid(lines = FALSE) +
    tm_xlab("Longitude") +
    tm_ylab("Latitude")

tmap_save(pep725_station_map, 
    filename = file.path("Output/Assets", "pep725_station_map.png"),
    width = 1600, height = 1800, dpi = 300)





# png(file.path("Output/Assets", "site_map.png"), 
#     width = 5200, height = 2000, res = 300)

# # Plot the main maps
# tmap_arrange(flux_sites_map, hb_map, pep725_station_map, ncol = 3)
# print(hf_hb_map, vp = viewport(0.6, 0.4, width = 0.3, height = 0.3))
# print(us_base_map, vp = viewport(0.6, 0.16, width = 0.2, height = 0.2))

# # Annotations
# grid.text("a", x = unit(0.05, "npc"), y = unit(0.98, "npc"))
# grid.text("b", x = unit(0.3, "npc"), y = unit(0.98, "npc"))
# grid.text("c", x = unit(0.75, "npc"), y = unit(0.98, "npc"))

# dev.off()

