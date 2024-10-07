# ******************************************************************************
# Summarize model fit comparison results from PEP725 + EOBS
# ******************************************************************************

species <- c("aesculus", "alnus", "betula", "fagus", "fraxinus", "quercus")

for (spec in species) {
    cat("Processing", spec, "\n")
    caches <- list.files("pipe/goodness-of-fit/pep", 
        pattern = spec,
        full.names = TRUE
    )

    pep_eobs <- lapply(caches, function(x) {
        mod <- readRDS(x)
        return(mod)
    })
    names(pep_eobs) <- lapply(caches, function(x) {
        strsplit(basename(x), "\\.")[[1]][1]
    })

    # out: model fit for PEP725 + E-OBS
    saveRDS(pep_eobs, file.path(
        "pipe/goodness-of-fit", 
        paste0("mod_", spec, "_pep725_eobs", ".Rds"))
    )
}
