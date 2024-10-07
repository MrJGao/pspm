# ******************************************************************************
# Process and visualize the result of the simulation.
# ******************************************************************************
source("src/base.R")
source("src/mod_hpc_pspm_sim_noise_base.R")
source("src/vis/vis_base.R")
library(data.table)
library(magrittr)



# Calculate the proportion of the true model identified as the best model.
CalBestModelProp <- function(dt, cur_model) {
    subdt <- dt[model == cur_model,]
    
    # For each noise level
    sum_dt <- by(subdt, subdt$noise, function (noise_dt) {
        sum_idx <- c()
        for (i in 1:nrow(noise_dt)) {
            irow <- noise_dt[i, ]
            idx <- which(irow[, .(tt_aicc, pa_aicc, sq_aicc, at_aicc, un_aicc)] - 
                irow$tt_aicc < 2
            )
            sum_idx <- c(sum_idx, idx)
        }
        return(data.table(
            noise = noise_dt$noise[1],
            tt_true = table(sum_idx)["1"] / length(sum_idx),
            pa_true = table(sum_idx)["2"] / length(sum_idx),
            sq_true = table(sum_idx)["3"] / length(sum_idx),
            at_true = table(sum_idx)["4"] / length(sum_idx),
            un_true = table(sum_idx)["5"] / length(sum_idx)
        ))
    }) %>%
        do.call(rbind, .)

    return(sum_dt)
}


DrawRect <- function(x, y, sum_dt) {
    size <- 0.8
    # Calculate TT proportions
    tt_xleft <- x - size / 2
    tt_xright <- (x - size / 2) + size * ifelse(
        is.na(sum_dt$tt_true), 0, sum_dt$tt_true
    )
    
    pa_xleft <- tt_xright
    pa_xright <- pa_xleft + size * ifelse(
        is.na(sum_dt$pa_true), 0, sum_dt$pa_true
    )
    
    sq_xleft <- pa_xright
    sq_xright <- sq_xleft + size * ifelse(
        is.na(sum_dt$sq_true), 0, sum_dt$sq_true
    )
    
    at_xleft <- sq_xright
    at_xright <- at_xleft + size * ifelse(
        is.na(sum_dt$at_true), 0, sum_dt$at_true
    )
    
    un_xleft <- at_xright
    un_xright <- un_xleft + size * ifelse(
        is.na(sum_dt$un_true), 0, sum_dt$un_true
    )

    # Draw rectangles
    # TT
    rect(xleft = tt_xleft, xright = tt_xright, 
        ybottom = y - size / 2, ytop = y + size / 2,
        col = jenna_pal[1]
    )
    # PA
    rect(xleft = pa_xleft, xright = pa_xright, 
        ybottom = y - size / 2, ytop = y + size / 2,
        col = jenna_pal[2]
    )
    # SQ
    rect(xleft = sq_xleft, xright = sq_xright, 
        ybottom = y - size / 2, ytop = y + size / 2,
        col = jenna_pal[3]
    )
    # AT
    rect(xleft = at_xleft, xright = at_xright, 
        ybottom = y - size / 2, ytop = y + size / 2,
        col = jenna_pal[4]
    )
    # UN
    rect(xleft = un_xleft, xright = un_xright, 
        ybottom = y - size / 2, ytop = y + size / 2,
        col = jenna_pal[5]
    )
    
}



sim_files <- list.files(file.path(hpc_dir, "Pipeline", "pspm_sim_noise"), 
    pattern = ".Rds$", full.names = TRUE
)

# ff <- sim_files[1]
dt <- lapply(sim_files, function(ff) {
    sim_fit <- readRDS(ff)

    # Retrieve simulation settings
    model <- sim_fit$sim$model
    noise <- sim_fit$sim$noise
    iter <- sim_fit$sim$iter

    # Calculate AICc
    tt_aicc <- AICc(sim_fit$TT$fitted_doy, sim_fit$obs, 3)$AICc
    pa_aicc <- AICc(sim_fit$PA$fitted_doy, sim_fit$obs, 9)$AICc
    sq_aicc <- AICc(sim_fit$SQ$fitted_doy, sim_fit$obs, 8)$AICc
    at_aicc <- AICc(sim_fit$AT$fitted_doy, sim_fit$obs, 5)$AICc
    un_aicc <- AICc(sim_fit$UN$fitted_doy, sim_fit$obs, 9)$AICc

    fit_dt <- data.table(
        model = model,
        noise = noise,
        iter = iter,
        tt_aicc = as.numeric(tt_aicc),
        pa_aicc = as.numeric(pa_aicc),
        sq_aicc = as.numeric(sq_aicc),
        at_aicc = as.numeric(at_aicc),
        un_aicc = as.numeric(un_aicc)
    )

    return(fit_dt)
})
dt <- do.call(rbind, dt)

{
    png("Output/pspm_sim_noise.png", width = 3000, height = 1500, res = 300)
    
    par(mar = c(5, 5, 3, 2))

    models <- c("TT", "PA", "SQ", "AT", "UN")
    plot(NA, xlim = c(-0.5, 10.5), ylim = c(0.5, 2.5),
        bty = "L", 
        xlab = "Variance of the observational noise (day)", ylab = "Model",
        yaxt = "n",
        cex.lab = 1.5, cex.axis = 1.2
    )
    axis(side = 2, las = 2, at = 1:2, labels = models[c(2, 4)], cex.lab = 1.2)

    # Draw perfect case
    # DrawRect(0, 1, data.table(tt_true = 1, pa_true = 0, sq_true = 0, 
    #     at_true = 0, un_true = 0)
    # )
    # DrawRect(0, 2, data.table(tt_true = 0, pa_true = 1, sq_true = 0, 
    #     at_true = 0, un_true = 0)
    # )
    # DrawRect(0, 3, data.table(tt_true = 0, pa_true = 0, sq_true = 1, 
    #     at_true = 0, un_true = 0)
    # )
    # DrawRect(0, 4, data.table(tt_true = 0, pa_true = 0, sq_true = 0, 
    #     at_true = 1, un_true = 0)
    # )
    
    # for (i in 1:length(models)) {
    #     if (models[i] %in% c("UN", "SQ", "TT")) {
    #         next
    #     }
    #     sum_dt <- CalBestModelProp(dt, cur_model = models[i])
    #     null <- lapply(sum_dt$noise, function(no) {
    #         DrawRect(no, i, sum_dt[noise == no])
    #     })
    # }

    DrawRect(0, 1, data.table(tt_true = 0, pa_true = 1, sq_true = 0, 
        at_true = 0, un_true = 0)
    )
    DrawRect(0, 2, data.table(tt_true = 0, pa_true = 0, sq_true = 0, 
        at_true = 1, un_true = 0)
    )
    sum_dt <- CalBestModelProp(dt, cur_model = "PA")
    null <- lapply(sum_dt$noise, function(no) {
        DrawRect(no, 1, sum_dt[noise == no])
    })
    sum_dt <- CalBestModelProp(dt, cur_model = "AT")
    null <- lapply(sum_dt$noise, function(no) {
        DrawRect(no, 2, sum_dt[noise == no])
    })


    legend(grconvertX(0.5, "ndc"), grconvertY(0.95, "ndc"), xjust = 0.5,
        ncol = 5, legend = models, fill = jenna_pal, xpd = NA, bty = "n",
        cex = 1.5
    )

    dev.off()
}



