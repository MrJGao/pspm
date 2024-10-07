#'******************************************************************************
#' Description: Compare models common functions.
#'******************************************************************************
library(data.table)
library(lubridate)


source("src/base.R")
null <- lapply(dir("src/hlp", full.names = TRUE), source)


FitUnModel <- function(data_li) {
    require(parallel)
    cl <- makeCluster(40)
    calls <- clusterCall(cl, function() {
        suppressWarnings({
            source("src/base.R")
            null <- lapply(dir("src/hlp", full.names = TRUE), source)
            source("src/mod_compare_base.R")
        })
    })
    clusterExport(cl, c("data_li"), envir = environment())


    mod_un_li <- clusterApplyLB(cl, 1:40, function(i) {
        set.seed(i)
        mod_un <- FitTheModel(
            par = NULL,
            model = UnifiedModel,
            cost_fun = CostRMSE,
            data = data_li,
            lower = c(1, 0, -10, -20, -100, 0, 0, -5, 10),
            upper = c(365, 10, 10, 5, 0, 20, 500, 0, 350),
            control = list(smooth = FALSE)
        )
        return(mod_un)
    })

    stopCluster(cl)

    best_idx <- sapply(mod_un_li, function(i) { i$optim$value })
    best_idx <- which.min(best_idx)
    mod_un <- mod_un_li[[best_idx]]

    return(mod_un)
}


FitUnModelSingle <- function(data_li) {
    require(parallel)
    
    mod_un_li <- lapply(1:30, function(i) {
        set.seed(i)
        mod_un <- FitTheModel(
            par = NULL,
            model = UnifiedModel,
            cost_fun = CostRMSE,
            data = data_li,
            lower = c(1, 0, -10, -20, -100, 0, 0, -5, 10),
            upper = c(365, 10, 10, 5, 0, 20, 500, 0, 350),
            control = list(smooth = FALSE)
        )
        return(mod_un)
    })

    best_idx <- sapply(mod_un_li, function(i) { i$optim$value })
    best_idx <- which.min(best_idx)
    mod_un <- mod_un_li[[best_idx]]

    return(mod_un)
}


#' Fit the classic spring phenology models and return the fit.
#' 
#' @param data_li The formatted data list.
#' @return Model fit variables.
FitCompareModels <- function(data_li, if_un_single = FALSE) {
    # TT
    mod_TT <- FitTheModel(
        model = ThermalTimeModel,
        cost_fun = CostRMSE,
        data = data_li
    )

    # PA
    mod_PA <- FitTheModel(
        model = ParallelModel,
        cost_fun = CostRMSE,
        data = data_li
    )

    # SQ
    mod_SQ <- FitTheModel(
        model = SequentialModel,
        cost_fun = CostRMSE,
        data = data_li
    )

    # AT
    mod_AT <- FitTheModel(
        model = AlternatingModel,
        cost_fun = CostRMSE,
        data = data_li
    )

    # UN
    # mod_UN <- FitTheModel(
    #     model = UnifiedModel,
    #     cost_fun = CostRMSE,
    #     data = data_li
    # )
     
    if (if_un_single == TRUE) {
        mod_UN <- FitUnModelSingle(data_li)
    } else {
        mod_UN <- FitUnModel(data_li)
    }

    return(list(
        obs = data_li$transition_dates,
        TT = mod_TT, PA = mod_PA, SQ = mod_SQ, AT = mod_AT, UN = mod_UN
    ))
}



#' Split the training and testing datasets by x-fold cv `idx`
#' The `idx` value varies from 1 to x, which is the number of folds.
#' 
#' @param data_li The formatted data list.
#' @param idx The test index of x-fold cv data.
#' @param folds Number of folds.
#' @return A list containing training and testing datasets.
SplitTrainTest <- function(data_li, idx, folds) {
    # `idx` specified fold goes to testing, others go to training
    fold_size <- floor(length(data_li$site) / folds)
    sy_idx_test <- seq((idx - 1) * fold_size + 1, idx * fold_size)
    sy_idx_train <- (1:length(data_li$site))[-sy_idx_test]

    testing_li <- GetDataLiIndex(data_li, sy_idx_test)
    training_li <- GetDataLiIndex(data_li, sy_idx_train)

    return(list(train = training_li, test = testing_li))
}


#' For leave-one-site-out cross validation, the data list will be split to 
#' train and test lists.
#' 
#' @param data_li The formatted data list.
#' @param siteID The site that will be left out.
#' @return A split list containing train and test lists.
LeaveOneSiteOut <- function(data_li, siteID) {
    sy_idx_test <- grep(siteID, data_li$site)
    sy_idx_train <- (1:length(data_li$site))[-sy_idx_test]

    training_li <- GetDataLiIndex(data_li, sy_idx_train)
    testing_li <- GetDataLiIndex(data_li, sy_idx_test)

    return(list(train = training_li, test = testing_li))
}

#' Cross-validation spring phenology models and return the fit.
#' To make it fair, the iteration times are all set to 1e5.
#' 
#' @param split_li The formatted data list that contains training and 
#' testing lists.
#' @return Model cv results, stored as anomalies.
CvCompareModels <- function(split_li, if_un_single = FALSE) {
    # TT
    mod_TT <- FitTheModel(
        model = ThermalTimeModel,
        cost_fun = CostRMSE,
        data = split_li$train
    )
    pred_TT <- do.call(ThermalTimeModel,
        list(par = mod_TT$optim$par, data = split_li$test))

    # PA
    mod_PA <- FitTheModel(
        model = ParallelModel,
        cost_fun = CostRMSE,
        data = split_li$train
    )
    pred_PA <- do.call(ParallelModel,
        list(par = mod_PA$optim$par, data = split_li$test))


    # SQ
    mod_SQ <- FitTheModel(
        model = SequentialModel,
        cost_fun = CostRMSE,
        data = split_li$train
    )
    pred_SQ <- do.call(SequentialModel,
        list(par = mod_SQ$optim$par, data = split_li$test))


    # AT
    mod_AT <- FitTheModel(
        model = AlternatingModel,
        cost_fun = CostRMSE,
        data = split_li$train
    )
    pred_AT <- do.call(AlternatingModel,
        list(par = mod_AT$optim$par, data = split_li$test))


    # UN
    # mod_UN <- FitTheModel(
    #     model = UnifiedModel,
    #     cost_fun = CostRMSE,
    #     data = split_li$train
    # )
    if (if_un_single == TRUE) {
        mod_UN <- FitUnModelSingle(split_li$train)
    } else {
        mod_UN <- FitUnModel(split_li$train)
    }

    pred_UN <- do.call(UnifiedModel,
        list(par = mod_UN$optim$par, data = split_li$test))


    return(list(
        obs = split_li$test$transition_dates,
        site = split_li$test$site,
        year = split_li$test$year,
        TT = pred_TT, PA = pred_PA, SQ = pred_SQ, AT = pred_AT, UN = pred_UN
    ))
}


#' Cross-validation spring phenology models and return the fit.
#' To make it fair, the iteration times are all set to 1e5.
#'
#' @param split_li The formatted data list that contains training and
#' testing lists.
#' @return Model cv results, stored as anomalies.
CvCompareModels2 <- function(split_li, ...) {
    # TT
    mod_TT <- FitTheModel(
        model = ThermalTimeModel,
        cost_fun = CostRMSE,
        data = split_li$train
    )
    pred_TT <- do.call(
        ThermalTimeModel,
        list(par = mod_TT$optim$par, data = split_li$test)
    )

    # PA
    mod_PA <- FitTheModel(
        model = ParallelModel,
        cost_fun = CostRMSE,
        data = split_li$train
    )
    pred_PA <- do.call(
        ParallelModel,
        list(par = mod_PA$optim$par, data = split_li$test)
    )

    # SQ
    mod_SQ <- FitTheModel(
        model = SequentialModel,
        cost_fun = CostRMSE,
        data = split_li$train
    )
    pred_SQ <- do.call(
        SequentialModel,
        list(par = mod_SQ$optim$par, data = split_li$test)
    )

    # AT
    mod_AT <- FitTheModel(
        model = AlternatingModel,
        cost_fun = CostRMSE,
        data = split_li$train
    )
    pred_AT <- do.call(
        AlternatingModel,
        list(par = mod_AT$optim$par, data = split_li$test)
    )

    # UN
    mod_UN <- FitTheModel(
        model = UnifiedModel,
        cost_fun = CostRMSE,
        data = split_li$train
    )
    pred_UN <- do.call(
        UnifiedModel,
        list(par = mod_UN$optim$par, data = split_li$test)
    )


    return(list(
        obs = split_li$test$transition_dates,
        site = split_li$test$site,
        year = split_li$test$year,
        TT = pred_TT, PA = pred_PA, SQ = pred_SQ, AT = pred_AT, UN = pred_UN
    ))
}


#' Cross-validation for a single model
#' 
#' @param data_li The formatted data list.
#' @param folds Number of folds.
#' @param model The spring phenology model to use.
#' @return A data table containing `siteID`, `PhenoYear`, `SOS`, and `pred_SOS`.
#' @export
CvModelFit <- function(data_li, folds, model, control = NULL) {
    cv_dt <- NULL
    for (i in 1:folds) {
        # Split train and test
        split_li <- SplitTrainTest(data_li, i, folds = folds)
        train_li <- split_li$train
        test_li <- split_li$test

        # Fit the model
        # Get par ranges from the json file
        par_range <- GetModelParRange(
            model_name = as.character(substitute(model))
        )
        par <- par_range$init
        lower <- par_range$lower
        upper <- par_range$upper

        est_par <- GenSA::GenSA(
            par = par,
            model = model,
            fn = CostRMSE,
            data = train_li,
            lower = lower,
            upper = upper,
            control = NULL
        )

        est <- do.call(model, list(par = est_par$par, data = train_li))
        mod_train <- list(optim = est_par, fitted_doy = est)

        

        # Test
        if (as.character(substitute(model)) == "LIN") {
            T_mean <- apply(train_li$Ti, 2, function(Ti) {
                mean(Ti[est_par$par[1]:est_par$par[2]])
            })
            mod_lm <- lm(train_li$transition_dates ~ T_mean)

            test_T_mean <- apply(test_li$Ti, 2, function(Ti) {
                mean(Ti[est_par$par[1]:est_par$par[2]])
            })
            mod_test <- predict(mod_lm, newdata = data.frame(T_mean = test_T_mean))
        } else {
            mod_test <- do.call(
                model,
                list(par = mod_train$optim$par, data = test_li)
            )
        }
        cv_dt <- rbind(cv_dt, data.table(
            siteID = test_li$siteID,
            PhenoYear = test_li$year,
            SOS = test_li$transition_dates,
            pred_SOS = mod_test
        ))
    }

    return(cv_dt)
}
