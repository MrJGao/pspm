#'******************************************************************************
#' Description: Functions for creating spring phenology model diagnosis figures. 
#' These figures can also be used to illustrate how the models work.
#'******************************************************************************
source("src/base.R")
source("src/vis/vis_base.R")

library(RColorBrewer)


mycolor <- brewer.pal(8, "Set2")
temp_color <- mycolor[1] # temperature color
sf_color <- mycolor[2] # forcing state color
gc_color <- mycolor[3] # growth competence color
sc_color <- mycolor[4] # chilling state color
rc_color <- mycolor[5] # chilling rate color


#' Internal use. Plot out temperature records.
#' @param data Needed data storing daily mean, min, and max temperature.
#' @return A temperture plot.
PlotTemperature <- function(data, idx) {
    plot(NA,
        xlim = range(data$doy), 
        ylim = range(data$T_max[, idx], data$T_min[, idx]),
        type = "l", xlab = "DOY", ylab = expression(
            "Tempeature(" * ~ degree * C * ")"),
        xaxt = "n"
    )
    grid()
    polygon(
        x = c(data$doy, rev(data$doy)), 
        y = c(data$T_max[, idx], rev(data$T_min[, idx])),
        col = Transparent(temp_color, 0.3), border = NA
    )
    lines(data$doy, data$Ti[, idx], lwd = 2, col = temp_color)
}


#' Internal use. Add observed and estimated phenological transition dates 
#' (in day-of-year format) onto the figure.
#' @param est_doy Estimated transition day of year.
#' @param obs_doy Observed transition day of year.
#' @return NULL
MarkPhenoDates <- function(est_doy, obs_doy) {
    # pheno dates ------------------
    abline(v = obs_doy)
    abline(v = est_doy, lty = 2)
    # annoate a text
    text(est_doy, grconvertY(0.7, "npc", "user"),
        pos = 2,
        labels = bquote(Est. - Obs. == .(est_doy - obs_doy) ~ days)
    )
    # add a legend
    legend("bottomright", lty = c(1, 2), legend = c("Obs.", "Est."), bty = "n")
}


#' Thermal Time model diagnosis figure
#' @param par Parameter vector.
#' @param data Needed data list.
#' @param idx The observation index for diagnosis.
#' @param est_doy The estimated transition day of year.
#' @param sf Forcing state.
#' @return A diagnosis figure.
#' @export
DiagThermalTimeFig <- function(par, data, idx, est_doy, sf) {
    t0 <- data$doy[round(par[1])]
    T_base <- par[2]
    F_crit <- par[3]

    # temperature -----------------------
    par(mfrow = c(2, 1), mar = c(0, 3, 1, 2))
    plot(NA, 
        xlim = range(data$doy), ylim = c(range(data$Ti[, idx])), xaxt = "n", 
        xlab = "", ylab = expression("Temperature (" * ~ degree * C * ")")
    )
    grid()
    lines(data$doy, data$Ti[, idx], col = temp_color, lwd = 2)
    # base temperature
    abline(h = T_base, lty = 3)
    text(-80, T_base, pos = 1, 
        labels = bquote(T[base] == .(round(T_base, 2)) ~ degree * C))
    
    # forcing state ------------------------
    par(mar = c(3, 3, 0, 2))
    plot(NA, xlim = range(data$doy), ylim = range(sf), xlab = "DOY", 
        ylab = "Forcing unit")
    grid()
    lines(data$doy, sf, lwd = 2, col = sf_color)
    # t0
    segments(t0, grconvertY(0, "npc", "user"), 
        t0, grconvertY(0.1, "npc", "user")
    )
    text(t0, grconvertY(0.1, "npc", "user"), 
        pos = 2, labels = bquote(t[0] == .(t0))
    )
    # forcing critical value
    abline(h = F_crit, lty = 3)
    text(-80, F_crit, pos = 3, labels = bquote(F[crit] == .(round(F_crit, 2))))
    
    # pheno dates ------------------
    MarkPhenoDates(est_doy, data$transition_dates[idx])
}


#' Parallel model diagnosis figure.
#' @param par Parameter vector.
#' @param data Needed data list.
#' @param idx The observation index for diagnosis.
#' @param est_doy The estimated transition day of year.
#' @param sf Forcing state.
#' @param sc Chilling state.
#' @param gc Growth competence.
#' @return A diagnosis figure.
#' @export
DiagParallelFig <- function(par, data, idx, est_doy, sf, sc, gc) {
    t0 <- data$doy[round(par[1])]
    t0_chill <- data$doy[round(par[2])]
    T_base <- par[3]
    T_opt <- par[4]
    T_min <- par[5]
    T_max <- par[6]
    C_min <- par[7]
    F_crit <- par[8]
    C_req <- par[9]


    # temperature ---------------------
    par(mfrow = c(4, 1), mar = c(0, 3, 1, 1), mgp = c(1.5, 0.5, 0))
    PlotTemperature(data, idx)
    # base temperature
    abline(h = T_base, lty = 3)
    text(-80, T_base, pos = 1, 
        labels = bquote(T[base] == .(round(T_base, 2)) ~ degree * C)
    )
    # min temperature
    segments(grconvertX(1, "npc", "user"), T_min, 
        grconvertX(0.95, "npc", "user"), T_min, lty = 3
    )
    text(grconvertX(0.95, "npc", "user"), T_min, pos = 2, 
        labels = bquote(T[min] == .(round(T_min, 2)) ~ degree * C))
    # max temperature
    segments(grconvertX(1, "npc", "user"), T_max, 
        grconvertX(0.95, "npc", "user"), T_max, lty = 3
    )
    text(grconvertX(0.95, "npc", "user"), T_max, pos = 2, 
        labels = bquote(T[max] == .(round(T_max, 2)) ~ degree * C))
    # optimal temperature
    segments(grconvertX(0, "npc", "user"), T_opt, 
        grconvertX(0.05, "npc", "user"), T_opt, lty = 3
    )
    text(grconvertX(0.05, "npc", "user"), T_opt, pos = 4, 
        labels = bquote(T[opt] == .(round(T_opt, 2)) ~ degree * C)
    )


    # growth competence function ------------------
    par(mar = c(0, 3, 0, 1))
    plot(NA, xlim = range(data$doy), ylim = c(0, 1), xaxt = "n", 
        ylab = "Growth competence")
    grid()
    lines(data$doy, gc, col = gc_color, lwd = 2)
    # minimum growth competence unchilled buds can have
    text(grconvertX(0.06, "npc", "user"), C_min, pos = 3, 
        labels = bquote(C[min] == .(round(C_min, 2))), col = gc_color)


    # chilling state ---------------------
    plot(NA, xlim = range(data$doy), ylim = range(sc, C_req), xaxt = "n", 
        ylab = "Chilling state")
    grid()
    lines(data$doy, sc, col = sc_color, lwd = 2)
    abline(h = C_req, lty = 3, col = sc_color)
    text(-80, C_req, pos = 1, labels = bquote(C[req] == .(round(C_req, 2))), 
        col = sc_color)


    # forcing unit ------------------
    par(mar = c(3, 3, 0, 1))
    plot(NA, xlim = range(data$doy), ylim = range(sf), xlab = "DOY", 
        ylab = "Forcing unit")
    grid()
    lines(data$doy, sf, lwd = 2, col = sf_color)
    # t0
    segments(t0, grconvertY(0, "npc", "user"), 
        t0, grconvertY(0.1, "npc", "user")
    )
    text(t0, grconvertY(0.1, "npc", "user"), 
        pos = 2, labels = bquote(t[0] == .(t0))
    )
    # t0_chill
    segments(t0_chill, grconvertY(0, "npc", "user"), t0_chill, 
        grconvertY(0.1, "npc", "user"))
    text(t0_chill, grconvertY(0.1, "npc", "user"), pos = 3, 
        labels = bquote(t[0][italic(chill)] == .(t0_chill)))
    # forcing critical value
    abline(h = F_crit, lty = 3)
    text(-80, F_crit, pos = 3, labels = bquote(F[crit] == .(round(F_crit, 2))))

    # pheno dates ------------------
    MarkPhenoDates(est_doy, data$transition_dates[idx])

}


#' Sequential model diagnosis figure.
#' @param par Parameter vector.
#' @param data Needed data list.
#' @param idx The observation index for diagnosis.
#' @param est_doy The estimated transition day of year.
#' @param sf Forcing state.
#' @param sc Chilling state.
#' @param gc Growth competence.
#' @return A diagnosis figure.
#' @export
DiagSequentialFig <- function(par, data, idx, est_doy, sf, sc, gc) {
    t0 <- data$doy[round(par[1])]
    t0_chill <- data$doy[round(par[2])]
    T_base <- par[3]
    T_opt <- par[4]
    T_min <- par[5]
    T_max <- par[6]
    F_crit <- par[7]
    C_req <- par[8]


    # plot temperature ---------------------
    par(mfrow = c(4, 1), mar = c(0, 3, 1, 1), mgp = c(1.5, 0.5, 0))
    PlotTemperature(data, idx)
    # base temperature
    abline(h = T_base, lty = 3)
    text(-80, T_base, pos = 1, 
        labels = bquote(T[base] == .(round(T_base, 2)) ~ degree * C)
    )
    # min temperature
    segments(grconvertX(1, "npc", "user"), T_min, 
        grconvertX(0.95, "npc", "user"), T_min, lty = 3
    )
    text(grconvertX(0.95, "npc", "user"), T_min, pos = 2, 
        labels = bquote(T[min] == .(round(T_min, 2)) ~ degree * C))
    # max temperature
    segments(grconvertX(1, "npc", "user"), T_max, 
        grconvertX(0.95, "npc", "user"), T_max, lty = 3
    )
    text(grconvertX(0.95, "npc", "user"), T_max, pos = 2, 
        labels = bquote(T[max] == .(round(T_max, 2)) ~ degree * C))
    # optimal temperature
    segments(grconvertX(0, "npc", "user"), T_opt, 
        grconvertX(0.05, "npc", "user"), T_opt, lty = 3
    )
    text(grconvertX(0.05, "npc", "user"), T_opt, pos = 4, 
        labels = bquote(T[opt] == .(round(T_opt, 2)) ~ degree * C)
    )


    # chilling state --------------------------------
    par(mar = c(0, 3, 0, 1))
    plot(NA, xlim = range(data$doy), ylim = range(sc, C_req), xaxt = "n", 
        ylab = "Chilling state")
    grid()
    lines(data$doy, sc, col = sc_color, lwd = 2)
    abline(h = C_req, lty = 3, col = sc_color)
    text(-80, C_req, pos = 1, labels = bquote(C[req] == .(round(C_req, 2))), 
        col = sc_color)
    # t0_chill
    segments(t0_chill, grconvertY(0, "npc", "user"), t0_chill, 
        grconvertY(0.1, "npc", "user"))
    text(t0_chill, grconvertY(0.1, "npc", "user"), pos = 3, 
        labels = bquote(t[0][~italic(chill)] == .(t0_chill)))


    # growth competence ------------------------------------
    plot(NA, xlim = range(data$doy), ylim = c(0, 1), ylab = "Growth competence")
    lines(data$doy, gc, col = gc_color, lwd = 2)


    # forcing state -------------------------------------
    par(mar = c(3, 3, 0, 1))
    plot(NA, xlim = range(data$doy), ylim = range(sf), xlab = "DOY", 
        ylab = "Forcing unit")
    grid()
    lines(data$doy, sf, lwd = 2, col = sf_color)
    # t0
    segments(t0, grconvertY(0, "npc", "user"), t0, 
        grconvertY(0.1, "npc", "user")
    )
    text(t0, grconvertY(0.1, "npc", "user"), 
        pos = 2, 
        labels = bquote(t[0] == .(t0))
    )
    # forcing critical value
    abline(h = F_crit, lty = 2)
    text(-80, F_crit, pos = 3, labels = bquote(F[crit] == .(round(F_crit, 2))))

    # pheno dates ------------------
    MarkPhenoDates(est_doy, data$transition_dates[idx])
}


#' Sequential model (J's version) diagnosis figure.
#' @param par Parameter vector.
#' @param data Needed data list.
#' @param idx The observation index for diagnosis.
#' @param est_doy The estimated transition day of year.
#' @param sf Forcing state.
#' @param sc Chilling state.
#' @param gc Growth competence.
#' @param t0_frc The DOY on which forcing accumulation starts.
#' @return A diagnosis figure.
#' @export
DiagSequentialFig_J <- function(par, data, idx, est_doy, sf, sc, gc, t0_frc) {
    t0 <- data$doy[round(par[1])]
    C_req <- par[2]
    T_base <- par[3]
    T_opt <- par[4]
    T_min <- par[5]
    T_max <- par[6]
    F_crit <- par[7]


    # plot temperature ---------------------
    par(mfrow = c(4, 1), mar = c(0, 3, 1, 1), mgp = c(1.5, 0.5, 0))
    PlotTemperature(data, idx)
    # base temperature
    abline(h = T_base, lty = 3)
    text(-80, T_base, pos = 1, 
        labels = bquote(T[base] == .(round(T_base, 2)) ~ degree * C)
    )
    # min temperature
    segments(grconvertX(1, "npc", "user"), T_min, 
        grconvertX(0.95, "npc", "user"), T_min, lty = 3
    )
    text(grconvertX(0.95, "npc", "user"), T_min, pos = 2, 
        labels = bquote(T[min] == .(round(T_min, 2)) ~ degree * C)
    )
    # max temperature
    segments(grconvertX(1, "npc", "user"), T_max, 
        grconvertX(0.95, "npc", "user"), T_max, lty = 3
    )
    text(grconvertX(0.95, "npc", "user"), T_max, pos = 2, 
        labels = bquote(T[max] == .(round(T_max, 2)) ~ degree * C)
    )
    # optimal temperature
    segments(grconvertX(0, "npc", "user"), T_opt, 
        grconvertX(0.05, "npc", "user"), T_opt, lty = 3
    )
    text(grconvertX(0.05, "npc", "user"), T_opt, pos = 4, 
        labels = bquote(T[opt] == .(round(T_opt, 2)) ~ degree * C)
    )


    # chilling state --------------------------------
    par(mar = c(0, 3, 0, 1))
    plot(NA, xlim = range(data$doy), ylim = range(sc, C_req), xaxt = "n", 
        ylab = "Chilling state")
    grid()
    lines(data$doy, sc, col = sc_color, lwd = 2)
    abline(h = C_req, lty = 3, col = sc_color)
    text(-80, C_req, pos = 3, labels = bquote(C[req] == .(round(C_req, 2))), 
        col = sc_color)
    # t0
    segments(t0, grconvertY(0, "npc", "user"), 
        t0, grconvertY(0.1, "npc", "user")
    )
    text(t0, grconvertY(0.1, "npc", "user"), 
        pos = 2, labels = bquote(t[0] == .(t0))
    )

    # growth competence ------------------------------------
    plot(NA, xlim = range(data$doy), ylim = c(0, 1), ylab = "Growth competence")
    grid()
    lines(data$doy, gc, col = gc_color, lwd = 2)


    # forcing state -------------------------------------
    par(mar = c(3, 3, 0, 1))
    plot(NA, xlim = range(data$doy), ylim = range(sf), xlab = "DOY", 
        ylab = "Forcing unit")
    grid()
    lines(data$doy, sf, lwd = 2, col = sf_color)
    # t0_frc
    segments(t0_frc, grconvertY(0, "npc", "user"), t0_frc, 
        grconvertY(0.1, "npc", "user"))
    text(t0_frc, grconvertY(0.1, "npc", "user"), pos = 2, 
        labels = bquote(t[0][~italic(forcing)] == .(t0_frc)))
    # forcing critical value
    abline(h = F_crit, lty = 2)
    text(-80, F_crit, pos = 3, labels = bquote(F[crit] == .(round(F_crit, 2))))

    # pheno dates ------------------
    MarkPhenoDates(est_doy, data$transition_dates[idx])
}


#' Alternating model diagnosis figure.
#' @param par Parameter vector.
#' @param data Needed data list.
#' @param idx The observation index for diagnosis.
#' @param est_doy The estimated transition day of year.
#' @return 
#' @export
DiagAlternatingFig <- function(par, data, idx, est_doy, sf, sfc, sc, rc) {
    # par = mod_AT$optim$par
    # data = site_li
    # est_doy = mod_AT$fitted_doy[idx]
    # sf = mod_AT_fit$Sf[, idx]
    # sfc = mod_AT_fit$Sfc[, idx]
    # sc = mod_AT_fit$Sc[, idx]
    # rc = mod_AT_fit$rc[, idx]

    # retrieve values from par
    t0 <- data$doy[round(par[1])]
    T_base <- par[2]
    a <- par[3]
    b <- par[4]
    c <- par[5]

    
    # temperature --------------------------------
    par(mfrow = c(4, 1), mar = c(0, 3, 1, 1), mgp = c(1.5, 0.5, 0))
    PlotTemperature(data, idx)
    # base temperature
    abline(h = T_base, lty = 3)
    text(-80, T_base, pos = 1, 
        labels = bquote(T[base] == .(round(T_base, 2)) ~ degree * C))


    # chilling rate -------------------------------
    par(mar = c(0, 3, 0, 1))
    plot(NA, xlim = range(data$doy), ylim = range(rc), xlab = "DOY", 
        ylab = "Chilling state")
    grid()
    lines(data$doy, rc, lwd = 2, col = sc_color)
    # t0
    segments(t0, grconvertY(0, "npc", "user"), 
        t0, grconvertY(0.1, "npc", "user")
    )
    text(t0, grconvertY(0.1, "npc", "user"), 
        pos = 2, labels = bquote(t[0] == .(t0))
    )


    # chilling state -------------------------------
    par(mar = c(0, 3, 0, 1))
    plot(NA, xlim = range(data$doy), ylim = range(sc), xlab = "DOY", 
        ylab = "Chilling state")
    grid()
    lines(data$doy, sc, lwd = 2, col = sc_color)
    # t0
    segments(t0, grconvertY(0, "npc", "user"), 
        t0, grconvertY(0.1, "npc", "user")
    )
    text(t0, grconvertY(0.1, "npc", "user"), 
        pos = 2, labels = bquote(t[0] == .(t0))
    )


    # forcing state -----------------------------------
    par(mar = c(3, 3, 0, 1))
    plot(NA, xlim = range(data$doy), ylim = range(sf, sfc), xlab = "DOY", 
        ylab = "Forcing state")
    grid()
    lines(data$doy, sf, lwd = 2, col = Transparent(sf_color, 0.3))
    lines(data$doy, sfc, lwd = 2, col = sf_color)
    # t0
    segments(t0, grconvertY(0, "npc", "user"), 
        t0, grconvertY(0.1, "npc", "user")
    )
    text(t0, grconvertY(0.1, "npc", "user"), 
        pos = 2, labels = bquote(t[0] == .(t0))
    )

    # pheno dates ----------------------
    MarkPhenoDates(est_doy, data$transition_dates[idx])
}