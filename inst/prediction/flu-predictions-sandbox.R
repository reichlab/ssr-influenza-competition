### post-process ssr fits to obtain predictions
library(lubridate)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)
library(ssrFlu)
library(cdcfluview)

locations <- c("ili_national", paste0("ili_region", 1:10))
pred_hzns <- 1:30

#######################################
## collect fits for a given location ##
#######################################

##for(location in locations) {
location <-  "ili_national" 
ssr_fits_by_prediction_horizon_limit <- lapply(pred_hzns,
                                               function(phl) {
                                                   file_name <- paste0(
                                                       "inst/estimation/2015-cdc-flu-competition/fit-competition-ssr-ph",
                                                       phl,
                                                       "-",
                                                       location,
                                                       ".Rdata")
                                                   read_env <- new.env()
                                                   load(file_name, envir=read_env)
                                                   return(read_env$ssr_fit)
                                               })
names(ssr_fits_by_prediction_horizon_limit) <- paste0("phl", pred_hzns)

## assemble data set

if(identical(location, "ili_national")) {
    usflu <- get_flu_data("national", "ilinet", years=1997:2015)
    data <- transmute(usflu,
                      region.type = REGION.TYPE,
                      region = REGION,
                      year = YEAR,
                      season_week = WEEK,
                      season = paste0(year, "/", year+1),
                      total_cases = as.numeric(X..WEIGHTED.ILI))
    
    ## subset is different than for estimation: adding 1 to index to account for rows removed for lag?
    data <- data[(262+1):nrow(data),]
    ## add time column
    data$time_ind <- seq_len(nrow(data))
} else {
    data <- NULL
}


##########################################
## make something resembling a forecast ##
##########################################

pred_horizon <- 1
ssr_fit <- ssr_fits_by_prediction_horizon_limit[[pred_horizon]]
max_lag <- max(unlist(ssr_fit$lags_hat))

## what are we predicting?
prediction_data_inds <- (nrow(data)-max_lag) : nrow(data)

## this adapted from challenge-predictions.R ~ lines 124-129
tmp <- ssr_predict(ssr_fit,
                   prediction_data=data[prediction_data_inds, , drop=FALSE],
                   leading_rows_to_drop=max(unlist(ssr_fit$lags_hat)),
                   prediction_horizon=pred_horizon,
                   normalize_weights=TRUE)
