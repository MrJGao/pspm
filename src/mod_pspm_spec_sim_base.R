#'******************************************************************************
#' Description: Base variables and functions for the species-level simulation
#' for exploring process-based spring phenology models.
#' 
#' Using HPC to faster the process.
#'******************************************************************************

source("src/base.R")
source("src/hlp/hlp_chilling_calc.R")
source("src/hlp/hlp_model_cost_functions.R")
source("src/hlp/hlp_spring_pheno_models.R")
source("src/hlp/hlp_fit_spm.R")
source("src/mod_compare_base.R")

library(data.table)
library(magrittr)

library(parallel)


# Use temperature data at site "1B" as an example b/c it has 27 site-years data
if (Sys.info()["sysname"] == "Darwin") { # Mac
    hfhb_dt <- fread(file.path(gdir, "Data/ARD/hb_hf_pheno_temp.csv"))
} else if (.Platform$OS.type == "windows") { # Windows
    hfhb_dt <- fread(file.path(gdir, "Data/ARD/hb_hf_pheno_temp.csv"))
} else if (.Platform$OS.type == "unix") { # Linux
    hfhb_dt <- fread(file.path(hpc_dir, "Data/ARD/hb_hf_pheno_temp.csv"))
} else {
    stop("Something is wrong...")
}


site1_li <- FormatDataForPhenoModel(hfhb_dt[siteID == "1B"])


# ~ Simulate data ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Species A - TT
a_tt_true <- c(
    t0 = 140,
    T_base = 5,
    F_crit = 300
)
a_sos <- do.call(ThermalTimeModel, list(
    par = a_tt_true, data = site1_li
))

# Species B - PA
b_pa_true <- c(
    t0 = 150,
    t0_chill = 90,
    T_base = 2,
    T_opt = 7,
    T_min = 0,
    T_max = 10,
    C_min = 0.1,
    F_crit = 200,
    C_req = 150
)
b_sos <- do.call(ParallelModel, list(
    par = b_pa_true, data = site1_li
))

# Species C - SQ
c_sq_true <- c(
    t0 = 100,
    t0_chill = 60,
    T_base = 5,
    T_opt = 2,
    T_min = -15,
    T_max = 20,
    F_crit = 400,
    C_req = 5
)
c_sos <- do.call(SequentialModel, list(
    par = c_sq_true, data = site1_li
))

# Species D - AT
d_at_true <- c(
    t0 = 150,
    T_base = 3,
    a = 50,
    b = 5000,
    c = -0.05
)
d_sos <- do.call(AlternatingModel, list(
    par = d_at_true, data = site1_li
))

# Species E - UN
e_un_true <- c(
    tc = 200,
    a_c = 0.1,
    b_c = 3,
    c_c = 2.6,
    b_f = -1,
    c_f = 10,
    w = 200,
    k = -0.015,
    C_req = 70
)
e_sos <- do.call(UnifiedModel, list(
    par = e_un_true, data = site1_li
))

# ~ Format data lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For species A
a_sim_li <- copy(site1_li)
a_sim_li$transition_dates <- a_sos

# For species B
b_sim_li <- copy(site1_li)
b_sim_li$transition_dates <- b_sos

# For species C
c_sim_li <- copy(site1_li)
c_sim_li$transition_dates <- c_sos

# For species D
d_sim_li <- copy(site1_li)
d_sim_li$transition_dates <- d_sos

# For species E
e_sim_li <- copy(site1_li)
e_sim_li$transition_dates <- e_sos


# Make all data into a single list for parallel processing
com_li <- list(
    a_sim_li, b_sim_li, c_sim_li, d_sim_li, e_sim_li
)
names(com_li) <- LETTERS[1:5]




# ~ Simulate species at different site locations ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# At this location, temperature is different, but model parameters are exactly
# the same with the first true model.

# Another site locaiton
site2_li <- FormatDataForPhenoModel(hfhb_dt[siteID == "5T"])

# Species A - TT
a2_sos <- do.call(ThermalTimeModel, list(
    par = a_tt_true, data = site2_li
))


# Species B - PA
b2_sos <- do.call(ParallelModel, list(
    par = b_pa_true, data = site2_li
))


# Species C - SQ
c2_sos <- do.call(SequentialModel, list(
    par = c_sq_true, data = site2_li
))


# Species D - AT
d2_sos <- do.call(AlternatingModel, list(
    par = d_at_true, data = site2_li
))


# Species E - UN
e2_sos <- do.call(UnifiedModel, list(
    par = e_un_true, data = site2_li
))


# ~ Format data lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For species A
a2_sim_li <- copy(site2_li)
a2_sim_li$transition_dates <- a2_sos

# For species B
b2_sim_li <- copy(site2_li)
b2_sim_li$transition_dates <- b2_sos

# For species C
c2_sim_li <- copy(site2_li)
c2_sim_li$transition_dates <- c2_sos

# For species D
d2_sim_li <- copy(site2_li)
d2_sim_li$transition_dates <- d2_sos

# For species E
e2_sim_li <- copy(site2_li)
e2_sim_li$transition_dates <- e2_sos


# Make all data into a single list for parallel processing
com2_li <- list(
    a2_sim_li, b2_sim_li, c2_sim_li, d2_sim_li, e2_sim_li
)
names(com2_li) <- LETTERS[1:5]





# ~ Simulate species at site 2 with slight changes in parameters ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# All change T_base

# Species A - TT
a3_tt_true <- a_tt_true
a3_tt_true["T_base"] <- 4
a3_sos <- do.call(ThermalTimeModel, list(
    par = a3_tt_true, data = site2_li
))


# Species B - PA
b3_pa_true <- b_pa_true
b3_pa_true["T_base"] <- 1
b3_sos <- do.call(ParallelModel, list(
    par = b3_pa_true, data = site2_li
))


# Species C - SQ
c3_sq_true <- c_sq_true
c3_sq_true["T_base"] <- 0
c3_sos <- do.call(SequentialModel, list(
    par = c3_sq_true, data = site2_li
))


# Species D - AT
d3_at_true <- d_at_true
d3_at_true["T_base"] <- 2
d3_sos <- do.call(AlternatingModel, list(
    par = d3_at_true, data = site2_li
))


# Species E - UN
e3_un_true <- e_un_true
e3_un_true["c_f"] <- 5
e3_sos <- do.call(UnifiedModel, list(
    par = e_un_true, data = site2_li
))


# ~ Format data lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For species A
a3_sim_li <- copy(site2_li)
a3_sim_li$transition_dates <- a3_sos

# For species B
b3_sim_li <- copy(site2_li)
b3_sim_li$transition_dates <- b3_sos

# For species C
c3_sim_li <- copy(site2_li)
c3_sim_li$transition_dates <- c3_sos

# For species D
d3_sim_li <- copy(site2_li)
d3_sim_li$transition_dates <- d3_sos

# For species E
e3_sim_li <- copy(site2_li)
e3_sim_li$transition_dates <- e3_sos


# Make all data into a single list for parallel processing
com3_li <- list(
    a3_sim_li, b3_sim_li, c3_sim_li, d3_sim_li, e3_sim_li
)
names(com3_li) <- LETTERS[1:5]




# ~ Aggregating ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Equally mixing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# At site 1, aggregate species equally
m1_sos <- (a_sos + b_sos + c_sos + d_sos + e_sos) / 5
m1_li <- copy(site1_li)
m1_li$transition_dates <- m1_sos

# ~ Chilling mixing only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# At site 1, aggregate chilling species only
m2_sos <- (b_sos + c_sos + d_sos + e_sos) / 4
m2_li <- copy(site1_li)
m2_li$transition_dates <- m2_sos


# ~ TT + PA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# At site 1, aggregate 70% species A with 30 % species B
m3_sos <- 0.7 * a_sos + 0.3 * b_sos
m3_li <- copy(site1_li)
m3_li$transition_dates <- m3_sos

# At site 1, aggregate 50% species A with 50% species B
m4_sos <- 0.5 * a_sos + 0.5 * b_sos
m4_li <- copy(site1_li)
m4_li$transition_dates <- m4_sos

# At site 1, aggregate 30% species A with 70% species B
m5_sos <- 0.3 * a_sos + 0.7 * b_sos
m5_li <- copy(site1_li)
m5_li$transition_dates <- m5_sos


# ~ TT + SQ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# At site 1, aggregate 70% species A with 30 % species C
m6_sos <- 0.7 * a_sos + 0.3 * c_sos
m6_li <- copy(site1_li)
m6_li$transition_dates <- m6_sos


# At site 1, aggregate 50% species A with 50% species C
m7_sos <- 0.5 * a_sos + 0.5 * c_sos
m7_li <- copy(site1_li)
m7_li$transition_dates <- m7_sos

# At site 1, aggregate 30% species A with 70% species C
m8_sos <- 0.3 * a_sos + 0.7 * c_sos
m8_li <- copy(site1_li)
m8_li$transition_dates <- m8_sos


# ~ TT + AT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# At site 1, aggregate 70% species A with 30 % species D
m9_sos <- 0.7 * a_sos + 0.3 * d_sos
m9_li <- copy(site1_li)
m9_li$transition_dates <- m9_sos

# At site 1, aggregate 50% species A with 50% species D
m10_sos <- 0.5 * a_sos + 0.5 * d_sos
m10_li <- copy(site1_li)
m10_li$transition_dates <- m10_sos

# At site 1, aggregate 30% species A with 70% species D
m11_sos <- 0.3 * a_sos + 0.7 * d_sos
m11_li <- copy(site1_li)
m11_li$transition_dates <- m11_sos


# ~ TT + UN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# At site 1, aggregate 70% species A with 30 % species E
m12_sos <- 0.7 * a_sos + 0.3 * e_sos
m12_li <- copy(site1_li)
m12_li$transition_dates <- m12_sos


# At site 1, aggregate 50% species A with 50% species E
m13_sos <- 0.5 * a_sos + 0.5 * e_sos
m13_li <- copy(site1_li)
m13_li$transition_dates <- m13_sos


# At site 1, aggregate 30% species A with 70% species E
m14_sos <- 0.3 * a_sos + 0.7 * e_sos
m14_li <- copy(site1_li)
m14_li$transition_dates <- m14_sos


mix_li <- list(
    m1_li, m2_li, m3_li, m4_li, m5_li,
    m6_li, m7_li, m8_li, m9_li, m10_li, m11_li, 
    m12_li, m13_li, m14_li
)
names(mix_li) <- paste0("m", 1:length(mix_li))



