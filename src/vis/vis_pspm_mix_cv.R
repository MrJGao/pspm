#'******************************************************************************
#' Description: Visualize the result of PSPM cross-validation simulated mixed
#' data.
#'******************************************************************************
library(data.table)
library(magrittr)
library(gridExtra)

source("src/base.R")

file_dir <- file.path(hpc_local_dir, "Pipeline", "pspm_mix_cv")



# ~ Goodness-of-fit ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gfit_rds <- list.files(file_dir,
    pattern = "pspm_spec_mix_fit_\\d+",
    full.names = TRUE
)

gfit_rmse_dt <- lapply(seq_along(gfit_rds), function(i) {
    i_gfit <- readRDS(gfit_rds[i])
    i_rmse <- lapply(seq_along(i_gfit), function(j) {
        x <- i_gfit[[j]]
        case <- j
        tt_rmse <- x$fit$TT$optim$value
        pa_rmse <- x$fit$PA$optim$value
        sq_rmse <- x$fit$SQ$optim$value
        at_rmse <- x$fit$AT$optim$value
        un_rmse <- x$fit$UN$optim$value

        return(data.table(
            iter = i,
            case = j,
            spec = i_gfit[[1]]$spec,
            tt_rmse, pa_rmse, sq_rmse, at_rmse, un_rmse
        ))
    }) %>%
        do.call(rbind, .)

    return(i_rmse)
}) %>%
    do.call(rbind, .)

# Calculate the base model for each iteration and each case
gfit_rmse_dt[,
    best := which.min(.SD[, .(
        tt_rmse, pa_rmse, sq_rmse,
        at_rmse, un_rmse
    )]),
    by = .(iter, case)
]


tab <- table(gfit_rmse_dt[, .(case, best)]) %>%
    data.table() %>%
    .[, lapply(.SD, as.numeric)]





# ~ Cross-validation ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ParseCvRds <- function(rds_files, spec_name) {
    res_dt <- lapply(seq_along(rds_files), function(i) {
        res <- readRDS(rds_files[i])
        res_dt <- lapply(res, function(yr) {
            fit <- data.frame(yr$fit)
            return(fit)
        }) %>%
            do.call(rbind, .) %>%
            data.table() %>%
            .[, ":="(spec = spec_name, iter = i)]

        # Remove NA
        res_dt[res_dt == 9999] <- NA

        return(res_dt)
    }) %>%
        do.call(rbind, .) %>%
        setcolorder(c("iter", "spec", "year"))

    # Calculate RMSE for each iteration for each model
    cv_rmse_dt <- res_dt[,
        .(
            tt_rmse = RMSE(.SD$TT, .SD$obs),
            pa_rmse = RMSE(.SD$PA, .SD$obs),
            sq_rmse = RMSE(.SD$SQ, .SD$obs),
            at_rmse = RMSE(.SD$AT, .SD$obs),
            un_rmse = RMSE(.SD$UN, .SD$obs)
        ),
        by = .(iter)
    ]

    # Calculate the best model
    cv_rmse_dt[,
        best := which.min(unlist(.SD[, .(
            tt_rmse, pa_rmse, sq_rmse,
            at_rmse, un_rmse
        )])),
        by = .(iter)
    ]

    return(cv_rmse_dt)
}


cv_tab_dt <- data.table()
for (i in 1:14) {
    case_rds <- list.files(file_dir,
        pattern = paste0("pspm_spec_mix_cv_\\d+_", i, ".Rds"),
        full.names = TRUE
    )

    case_dt <- ParseCvRds(case_rds, paste0("m", i))


    case_tab <- table(case_dt$best) / nrow(case_dt) * 100

    case_dt <- rep(0, 5)
    case_dt[names(case_tab) %>% as.numeric()] <- case_tab
    case_dt <- data.table(case = i, t(case_dt))
    names(case_dt) <- c("case", "tt", "pa", "sq", "at", "un")
    
    cv_tab_dt <- rbind(cv_tab_dt, case_dt)
}



# ~ Make the figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
case_names <- c(
    "Equally mix", "Chilling mix",
    "70%TT+30%PA", "50%TT+50%PA", "30%TT+70%PA",
    "70%TT+30%SQ", "50%TT+50%SQ", "30%TT+70%SQ",
    "70%TT+30%AT", "50%TT+50%AT", "30%TT+70%AT",
    "70%TT+30%UN", "50%TT+50%UN", "30%TT+70%UN"
)

{

png("Output/Assets/sim_mix.png", width = 2000, height = 2600, res = 150)

# png("Output/Assets/sim_mix_dark.png", width = 3700, height = 4000, res = 300)
# par(bg = NA, fg = "white", col.axis = "white", col.lab = "white",
#     col.main = "white", col.sub = "white")

par(mfrow = c(5, 3), cex.lab = 1.5, cex.axis = 1.5, oma = c(0, 3, 0, 0))
layout(matrix(c(
    1, 2, 3,
    1, 2, 3,
    4:15
), nrow = 6, byrow = TRUE))
# Goodness-of-fit
par(mar = c(4, 11, 1, 3))
plot(NA,
    xlim = c(0.5, 5), ylim = c(1, 14),
    bty = "L",
    xlab = "", xaxt = "n",
    ylab = "", yaxt = "n"
)
grid()
axis(side = 1, at = 1:5, labels = c("TT", "PA", "SQ", "AT", "UN"))
axis(side = 2, at = 1:14, labels = case_names,
    las = 1
)
mtext(side = 2, text = "Case", line = 11)
points(tab[, .(best, case)], pch = 15, cex = 100 / tab$N * 1)
text(tab[N !=0 , .(best, case)], labels = tab[N !=0, paste0(N, "%")], 
    cex = 1.5, pos = 2, xpd = NA
)

# Cross-validation
par(mar = c(2, 3.5, 3, 3))
for (i in 1:11) {
    barplot(cv_tab_dt[case == i, 2:6] %>% unlist(),
        names.arg = c("TT", "PA", "SQ", "AT", "UN"),
        ylab = "Percentage (%)",
        xlab = "",
        xaxt = "n",
        xpd = NA
    )
    title(case_names[i], cex.main = 2)
}
for (i in 12:14) {
    barplot(cv_tab_dt[case == i, 2:6] %>% unlist(),
        names.arg = c("TT", "PA", "SQ", "AT", "UN"),
        ylab = "Percentage (%)",
        xlab = "",
        xpd = NA
    )
    title(case_names[i], cex.main = 2)
}

# Annotation
text(grconvertX(rep(c(0.02, 0.35, 0.67), 2), "ndc"), 
    grconvertY(rep(c(0.98, 0.79, 0.59, 0.39, 0.19), each = 3), "ndc"),
    labels = letters[1:15],
    xpd = NA,
    cex = 3
)

dev.off()

}


