#'******************************************************************************
#' Description: Download E-OBS daily mean temperature dataset.
#'******************************************************************************
library(httr)
library(rvest)
library(xml2)

source("src/base.R")

# in:
# The E-OBS data website url
url <- "https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php"

# EOBS/
out_dir <- "data/raw/EOBS"


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

the_page <- read_html(url)
# Find the table
the_table <- html_node(the_page, "#block-system-main") %>%
    html_nodes("table") %>% 
    .[[2]]
# Find the second row, which is version 24.0e
the_v240e_row <- html_nodes(the_table, "tr") %>% .[[2]]

# Find the link with `TG` which means `mean temperature`
the_ds_cell <- html_nodes(the_v240e_row, "td") %>% .[[2]]
items <- html_children(the_ds_cell) %>% html_text()
links <- html_children(the_ds_cell) %>% html_attr("href")
data_link_df <- data.frame(
    DataName = items, 
    Link = links, 
    stringsAsFactors = FALSE
)

tg_link <- data_link_df[data_link_df$DataName == "TG", "Link"]

# Find the elevation cell
the_elev_cell <- html_nodes(the_v240e_row, "td") %>% 
    .[[4]] %>%
    html_children()
elev_link <- html_attr(the_elev_cell, "href")


# out: e-obs tg
tg_file <- file.path(out_dir, basename(tg_link))

# Dl E-OBS TG data
response <- GET(tg_link, write_disk(
        path = tg_file,
        overwrite = TRUE
), progress())


# out: e-obs elev
elev_file <- file.path(out_dir, basename(elev_link))

# Dl E-OBS elevation
response <- GET(elev_link, write_disk(
    path = elev_file,
    overwrite = TRUE
), progress())

# NOTE: The E-OBS elevation link is WRONG by the time I download data! The 
# downloaded version 24.0e of elevation data file is actually a day of TG 
# temperature data. So, I downloaded version 23.1e elevation instead. 
# The file link should be fixed in the future.
