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
    
	
	## add log column
	data$log_total_cases <- log(data$total_cases + 1)
	
	## add week_start_date
	## note we changed reference to data$week in estimation file to data$season_week here
	char_dates <- paste(data$year, data$season_week, "1")
	data$week_start_date <- as.Date(char_dates, format="%Y %W %w")
	
	## remove week 53s
	data <- data[-which(is.na(data$week_start_date)),]
	## subset is different than for estimation: adding 1 to index to account for rows removed for lag?
	data <- data[(262+1):nrow(data),]
	
	## add smooth log column -- or not, since it is already pretty smooth...
	## (but this was done in estimation, so we need to do it here as well.)
	sm <- loess(log_total_cases ~ as.numeric(week_start_date), data=data, span= 26 / nrow(data))
	data$smooth_log_cases <- sm$fitted
	
    ## add time column
    data$time_ind <- seq_len(nrow(data))
} else {
    data <- NULL
}


##########################################
## make something resembling a forecast ##
##########################################

nsim <- 1000
pred_horizons <- 1:30
#pred_horizons <- 1:10  THIS ONE WORKS.
preds <- matrix(nrow=nsim*length(pred_horizons), ncol=3)
colnames(preds) <- c("hzn", "sim", "ili")
preds[,"hzn"] <- rep(pred_horizons, each=nsim)
preds[,"ili"] <- rep(1:nsim, times=length(pred_horizons))

for(i in 1:length(pred_horizons)){
    message(paste("horizon", i))
    pred_horizon <- pred_horizons[i]
    ssr_fit <- ssr_fits_by_prediction_horizon_limit[[pred_horizon]]
    max_lag <- max(unlist(ssr_fit$lags_hat))
    
    ## what are we predicting?
    prediction_data_inds <- (nrow(data)-max_lag) : nrow(data)
    
    ## update theta_est in ssr_fit object to also contain parameter values that were held fixed
    ## in this case, this includes only the period of the periodic kernel function
    if("time_ind_lag0" %in% names(ssr_fit$theta_hat)) {
        ssr_fit$theta_hat$time_ind_lag0 <- c(ssr_fit$theta_hat$time_ind_lag0,
                                             ssr_fit$ssr_control$theta_fixed$time_ind)
    }
    
    ## this adapted from challenge-predictions.R ~ lines 124-129
    tmp <- ssr_predict(ssr_fit,
                       prediction_data=data[prediction_data_inds, , drop=FALSE],
                       leading_rows_to_drop=max(unlist(ssr_fit$lags_hat)),
                       prediction_horizon=pred_horizon,
                       normalize_weights=TRUE)
    
    idx <- which(preds[,"hzn"]==i)
    preds[idx,"ili"] <- simulate_from_weighted_kde(1000, tmp)
}

preds_sum <- tbl_df(data.frame(preds)) %>%
    group_by(hzn) %>%
    summarize(median_ili = median(ili),
              p05 = quantile(ili, .05),
              p95 = quantile(ili, .95))


ggplot(preds_sum, aes(x=hzn)) + 
    geom_line(aes(y=median_ili)) +
    geom_ribbon(aes(ymin=p05, ymax=p95), alpha=.2)




## draw simulated values from the distribution encoded as a weighted kernel density estimate in the tmp object
## this uses the simulate_from_weighted_kde function above, which makes use of the
##  - weights: vector of the weight assigned to each kernel center
##  - centers: vector of centers for the kernels
## the simulate_from_weighted_kde function uses R's built in "density" function to estimate the predictive distribution bandwidth
# sample <- simulate_from_weighted_kde(1000, tmp)
# 
# par(mfrow = c(2, 1))
# plot(sample)
# hist(sample)

