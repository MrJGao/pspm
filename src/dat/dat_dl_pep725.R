
library(httr)
library(xml2)
library(rvest)
library(data.table)
source("src/base.R")


# in:
# PEP725 data website
pep_dt_url <- "http://www.pep725.eu/data_download/data_selection.php"
# PEP725 login url
pep_login_url <- "http://www.pep725.eu/login.php"
# PEP725 data selection url
pep_data_sel_url <- "http://www.pep725.eu/data_download/data_selection.php"


# out: data/raw/PEP725/
out_dir <- "data/raw/PEP725"


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data_selection <- httr::GET(pep_dt_url)
sel_html <- read_html(data_selection)
sel_node <- html_nodes(sel_html, "form select")
sel_items <- html_children(sel_node)

item_names <- html_text(sel_items)
item_numbers <- as.numeric(html_attr(sel_items, "value"))

# All species and their code
sel_options <- data.table(spec = item_names, code = item_numbers)


# Get the species Zhang et al., 2022 used.
# I can only make it work by changing some of the names.
# sel_options[spec %in% c("Aesculus hippocastanum", "Alnus glutinosa", 
# "Betula pendula", "Fagus sylvatica", "Fraxinus excelsior", "Quercus robur"), ]
target_species <- sel_options[spec %in% c(
    "Aesculus hippocastanum", "Alnus",
    "Betula", "Fagus", "Fraxinus excelsior", "Quercus robur (Q.peduncula)"
), ]

# Download data
account <- ParsePrivate("pep725")

# Login to the website
login <- POST(pep_login_url, body = list(
    email = account$email,
    pwd = account$password, submit = "Login"
), encode = "form")

# Download species data
species_html <- lapply(target_species$code, function(x) {
    POST(pep_data_sel_url,
        body = list(
            plant = x,
            submit1 = "Submit"
        ),
        encode = "form"
    )
})

# extract the links to download
species_links <- lapply(species_html, function(x) {
    read_html(x) %>%
        rvest::html_nodes("td a") %>%
        rvest::html_attr("href")
})
species_links <- do.call(c, species_links)

# Process downloading
lapply(species_links, function(link) {
    GET(link, write_disk(
        path = file.path(out_dir, sprintf(
                    "PEP725_%s.tar.gz",
                    strsplit(link, "=")[[1]][2]
                )),
        overwrite = TRUE
    ))
})
