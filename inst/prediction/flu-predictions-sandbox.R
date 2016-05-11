### post-process ssr fits to obtain predictions
library(ggplot2)
#library(ssrFlu)
library(cdcfluview)

library(lubridate)
library(reshape)
library(plyr)
library(dplyr)
library(tidyr)

library(copula)
library(mvtnorm)

library(kcde)

## Function borrowed from copula package and modified
## Find a matrix near mat that is positive definite
makePosDef <- function (mat, delta = 0.001) 
{
    while(min(eigen(mat)$values) < 10^{-6}) {
        decomp <- eigen(mat)
        Lambda <- decomp$values
        Lambda[Lambda < 0] <- delta
        Gamma <- decomp$vectors
        newmat <- Gamma %*% diag(Lambda) %*% t(Gamma)
        D <- 1/sqrt(diag(newmat))
        mat <- diag(D) %*% newmat %*% diag(D)
    }
    return(mat)
}


## set up -- these parameters determine the data set and which KCDE fit we will use
data_set <- "ili_national"
max_lag <- 1L
max_seasonal_lag <- 0L
filtering <- FALSE
differencing <- FALSE
seasonality <- TRUE
#bw_parameterization <- "diagonal"
bw_parameterization <- "full"

n_sims <- 10000

copula_save_path <- file.path("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
    data_set,
    "copula-estimation-results")
estimation_save_path <- file.path("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
    data_set,
    "estimation-results")
prediction_save_path <- file.path("/media/evan/data/Reich/infectious-disease-prediction-with-kcde/inst/results",
    data_set,
    "prediction-results")

case_descriptor <- paste0(
    data_set,
    "-max_lag_", max_lag,
    "-max_seasonal_lag_", max_seasonal_lag,
    "-filtering_", filtering,
    "-differencing_", differencing,
    "-seasonality_", seasonality,
    "-bw_parameterization_", bw_parameterization
)
file_name <- paste0("kcde-copula-fits-",
    case_descriptor,
    ".rds")
copula_fits <- readRDS(file = file.path(copula_save_path, file_name))
analysis_time_season_week_by_copula_fit <- unlist(lapply(copula_fits,
        function(copula_fit) {copula_fit$analysis_time_season_week}))


kcde_fits_by_prediction_horizon <- lapply(seq_len(52),
    function(prediction_horizon) {
        case_descriptor <- paste0(
            data_set,
            "-prediction_horizon_", prediction_horizon,
            "-max_lag_", max_lag,
            "-max_seasonal_lag_", max_seasonal_lag,
            "-filtering_", filtering,
            "-differencing_", differencing,
            "-seasonality_", seasonality,
            "-bw_parameterization_", bw_parameterization
        )
        readRDS(file.path(estimation_save_path,
                paste0("kcde_fit-", case_descriptor, ".rds")))
    })


if(identical(data_set, "ili_national")) {
    ## Load data for nationally reported influenza like illness
    library(cdcfluview)
    
    usflu <- get_flu_data("national", "ilinet", years=1997:2015)
    data <- transmute(usflu,
        region.type = REGION.TYPE,
        region = REGION,
        year = YEAR,
        week = WEEK,
        weighted_ili = as.numeric(X..WEIGHTED.ILI))
    
    ## Add time column.  This is used for calculating times to drop in cross-validation
    data$time <- ymd(paste(data$year, "01", "01", sep = "-"))
    week(data$time) <- data$week
    
    ## Add time_index column.  This is used for calculating the periodic kernel.
    ## Here, this is calculated as the number of days since some origin date (1970-1-1 in this case).
    ## The origin is arbitrary.
    data$time_index <- as.integer(data$time -  ymd(paste("1970", "01", "01", sep = "-")))
    
    ## Season column: for example, weeks of 2010 up through and including week 30 get season 2009/2010;
    ## weeks after week 30 get season 2010/2011
    data$season <- ifelse(
        data$week <= 30,
        paste0(data$year - 1, "/", data$year),
        paste0(data$year, "/", data$year + 1)
    )
    
    ## Season week column: week number within season
    data$season_week <- sapply(seq_len(nrow(data)), function(row_ind) {
            sum(data$season == data$season[row_ind] & data$time_index <= data$time_index[row_ind])
        })
    
    if(differencing) {
        data$weighted_ili_ratio <- data$weighted_ili / lag(data$weighted_ili, 52) 
        prediction_target_var <- "weighted_ili_ratio"
        train_seasons <- paste0(seq(from = 1998, to = 2009), "/", seq(from = 1999, to = 2010))
    } else {
        prediction_target_var <- "weighted_ili"
        train_seasons <- paste0(seq(from = 1997, to = 2009), "/", seq(from = 1998, to = 2010))
    }
    
    kernel_fn <- log_pdtmvn_kernel
    rkernel_fn <- rlog_pdtmvn_kernel
    
    variable_selection_method <- "all_included"
    crossval_buffer <- ymd("2010-01-01") - ymd("2009-01-01")
    
    season_length <- 33L
    analysis_seasons <- "2015/2016"
    first_analysis_time_season_week <- 10 # == week 40 of year
    last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
}


last_obs_week <- 17
last_obs_year <- 2016
analysis_time_ind <- which(data$year == last_obs_year & data$week == last_obs_week)
analysis_time_season <- data$season[analysis_time_ind]
analysis_time_season_week <- data$season_week[analysis_time_ind]

### simulate from copula that ties the marginal predictive distributions together

## get the right copula for analysis_time_season_week
predictive_copula_ind <- which(analysis_time_season_week_by_copula_fit == analysis_time_season_week)
copula_fit <- copula_fits[[predictive_copula_ind]]$copula_fit
predictive_copula <- copula_fit@copula

## simulate n_sims sequences from copula
max_prediction_horizon <-
    last_analysis_time_season_week + 1 -
    analysis_time_season_week
if(max_prediction_horizon < 4L) {
    ## get the right copula for analysis_time_season_week
    predictive_copula_ind <- which(analysis_time_season_week_by_copula_fit == analysis_time_season_week) -
        (max_prediction_horizon - 4L)
    copula_fit <- copula_fits[[predictive_copula_ind]]$copula_fit
    predictive_copula <- copula_fit@copula
    
    max_prediction_horizon <- 4L
}
sim_sequences <- matrix(NA, nrow = n_sims, ncol = max_prediction_horizon)
na_rows <- seq_len(n_sims)
#    while(length(na_rows) > 0) {
for(sim_ind in na_rows) {
    ## generate random parameters from estimation distribution
    ## Note that the variance estimate for parameters is low; at least we're
    ## accounting for some uncertainty though...
#    predictive_copula@parameters <- rmvnorm(1, copula_fit@estimate, sigma = copula_fit@var.est)[1, ]
#    if(identical(class(predictive_copula)[1], "normalCopula")) {
#        ## randomly generated parameters above may not yield a positive definite correlation matrix; correct
#        predictive_copula@parameters <- makePosDef(getSigma(predictive_copula))[1, 1 + seq_along(predictive_copula@parameters)]
#    }
#    
    ## simulate sequence from copula
    sim_sequences[sim_ind, ] <- rCopula(1, predictive_copula)[1, ]
}

## NA rows result when random parameters are problematic
#        na_rows <- which(apply(sim_sequences, 1, function(ss_row) any(is.na(ss_row))))
#    }

### get quantiles from marginal predictive distributions corresponding to
### values simulated from copula
#analysis_time_ind <- which(data$season == analysis_time_season &
#        data$season_week == analysis_time_season_week)
trajectory_samples <- matrix(NA, nrow = n_sims, ncol = max_prediction_horizon)
for(prediction_horizon in seq_len(max_prediction_horizon)) {
    trajectory_samples[, prediction_horizon] <-
        kcde_predict(
            p = sim_sequences[, prediction_horizon],
            n = 100000,
            kcde_fit = kcde_fits_by_prediction_horizon[[prediction_horizon]],
            prediction_data =
                data[seq_len(analysis_time_ind), , drop = FALSE],
            leading_rows_to_drop = 0L,
            trailing_rows_to_drop = 0L,
            additional_training_rows_to_drop = NULL,
            prediction_type = "quantile",
            log = TRUE
        )
    if(differencing) {
        trajectory_samples[, prediction_horizon] <-
            trajectory_samples[, prediction_horizon] *
            data[analysis_time_ind + prediction_horizon - 52, prediction_target_var]
    }
}



#locations <- c("ili_national", paste0("ili_region", 1:10))
#pred_hzns <- 1:30
#
########################################
### collect fits for a given location ##
########################################
#
###for(location in locations) {
#location <-  "ili_national" 
#ssr_fits_by_prediction_horizon_limit <- lapply(pred_hzns,
#                                               function(phl) {
#                                                   file_name <- paste0(
#                                                       "inst/estimation/2015-cdc-flu-competition/fit-competition-ssr-ph",
#                                                       phl,
#                                                       "-",
#                                                       location,
#                                                       ".Rdata")
#                                                   read_env <- new.env()
#                                                   load(file_name, envir=read_env)
#                                                   return(read_env$ssr_fit)
#                                               })
#names(ssr_fits_by_prediction_horizon_limit) <- paste0("phl", pred_hzns)
#
### assemble data set
#
#if(identical(location, "ili_national")) {
#    usflu <- get_flu_data("national", "ilinet", years=1997:2015)
#    data <- transmute(usflu,
#                      region.type = REGION.TYPE,
#                      region = REGION,
#                      year = YEAR,
#                      season_week = WEEK,
#                      season = paste0(year, "/", year+1),
#                      total_cases = as.numeric(X..WEIGHTED.ILI))
#    
#	
#	## add log column
#	data$log_total_cases <- log(data$total_cases + 1)
#	
#	## add week_start_date
#	## note we changed reference to data$week in estimation file to data$season_week here
#	char_dates <- paste(data$year, data$season_week, "1")
#	data$week_start_date <- as.Date(char_dates, format="%Y %W %w")
#	
#	## remove week 53s
#	data <- data[-which(is.na(data$week_start_date)),]
#	## subset is different than for estimation: adding 1 to index to account for rows removed for lag?
#	data <- data[(262+1):nrow(data),]
#	
#	## add smooth log column -- or not, since it is already pretty smooth...
#	## (but this was done in estimation, so we need to do it here as well.)
#	sm <- loess(log_total_cases ~ as.numeric(week_start_date), data=data, span= 26 / nrow(data))
#	data$smooth_log_cases <- sm$fitted
#	
#    ## add time column
#    data$time_ind <- seq_len(nrow(data))
#} else {
#    data <- NULL
#}


##########################################
## make something resembling a forecast ##
##########################################

filedate <- '20160506'
last_obs_week <- 17
last_obs_year <- 2016
nsim <- 10000
pred_horizons <- seq_len(ncol(trajectory_samples))
#pred_horizons <- 1:30
#pred_horizons <- 1:10  THIS ONE WORKS.
preds <- matrix(nrow=nsim*length(pred_horizons), ncol=3)
colnames(preds) <- c("hzn", "sim", "ili")
preds[,"hzn"] <- rep(pred_horizons, each=nsim)
preds[,"sim"] <- rep(1:nsim, times=length(pred_horizons))

for(i in 1:length(pred_horizons)){
#    message(paste("horizon", i))
#    pred_horizon <- pred_horizons[i]
#    ssr_fit <- ssr_fits_by_prediction_horizon_limit[[pred_horizon]]
#    max_lag <- max(unlist(ssr_fit$lags_hat))
#    
#    ## what are we predicting?
#    prediction_data_inds <- (nrow(data)-max_lag) : nrow(data)
#    
#    ## update theta_est in ssr_fit object to also contain parameter values that were held fixed
#    ## in this case, this includes only the period of the periodic kernel function
#    if("time_ind_lag0" %in% names(ssr_fit$theta_hat)) {
#        ssr_fit$theta_hat$time_ind_lag0 <- c(ssr_fit$theta_hat$time_ind_lag0,
#                                             ssr_fit$ssr_control$theta_fixed$time_ind)
#    }
#    
#    ## this adapted from challenge-predictions.R ~ lines 124-129
#    tmp <- ssr_predict(ssr_fit,
#                       prediction_data=data[prediction_data_inds, , drop=FALSE],
#                       leading_rows_to_drop=max(unlist(ssr_fit$lags_hat)),
#                       prediction_horizon=pred_horizon,
#                       normalize_weights=TRUE)
    
    idx <- which(preds[,"hzn"]==i)
    preds[idx,"ili"] <- trajectory_samples[, i]
#    preds[idx,"ili"] <- simulate_from_weighted_kde(nsim, tmp)
}

#upto_current_data_idx <- which(data$season=="2015/2016" & data$season_week==40):nrow(data)
#data_sims <- data_frame(week = rep(data[upto_current_data_idx, "season_week"], each=nsim), 
#    sim = rep(1:nsim, times=length(upto_current_data_idx)),
#    season_week = 201500 + 100*(week<30) + week,
#    ili = rep(data[upto_current_data_idx, "total_cases"], each=nsim))

upto_current_data_idx <- which(data$season=="2015/2016" & data$week==40):nrow(data)
data_sims <- data_frame(week = rep(data[upto_current_data_idx, "week"], each=nsim), 
                        sim = rep(1:nsim, times=length(upto_current_data_idx)),
                        season_week = 201500 + 100*(week<30) + week,
                        ili = rep(data[upto_current_data_idx, "weighted_ili"], each=nsim))

preds_df <- tbl_df(data.frame(preds)) %>%
    mutate(week = (last_obs_week + hzn - 1)%%52 + 1,
           season_week = 201500 + 100*(week<30) + week) %>%
    full_join(data_sims) %>%
    mutate(week_date = as.Date(paste0(season_week,00), format="%Y%W%w")) %>%
    arrange(season_week)

#########################################
## calculuate peak week probabilities  ##
#########################################

template_table <- data_frame(week=c(40:52, 1:20), season_week=c(201540:201552, 201601:201620))

peak_week <- preds_df %>% group_by(sim) %>%
    summarize(week = week[which.max(ili)]) %>%
    ungroup() %>% group_by(week) %>%
    summarize(peak_wk_totals = n()) %>%
    ungroup() %>% right_join(template_table) 
## replace zeroes
peak_week[is.na(peak_week)] <- 0
peak_week <- peak_week %>%
    mutate(peak_wk_totals = (peak_wk_totals+1),
           peak_wk_prob = peak_wk_totals/sum(peak_wk_totals))

write.csv(peak_week, file=paste0('inst/submissions/', filedate, '-peak-week.csv'))
    
############################################
## calculuate season onset probabilities  ##
############################################

## need to add no onset possibility?

template_onsets <- data_frame(first_week_year = c(rep(2015, 13), rep(2016, 20)),
                              first_week_num = c(40:52, 1:20))

seasonal_baseline <- 2.1
onsets <- preds_df %>% group_by(sim) %>%
    mutate(ili_lag1 = lag(ili, 1),
           ili_lag2 = lag(ili, 2),
           ili_lag3 = lag(ili, 3),
           onset = ili_lag1 > seasonal_baseline & ili_lag2 > seasonal_baseline & ili_lag3 > seasonal_baseline) %>%
    filter(onset) %>% 
    summarize(first_week = first(week_date, order_by=onset)-weeks(2)) %>%
    ungroup() %>%
    mutate(first_week_num = as.numeric(format(first_week, "%W")),
           first_week_year = as.numeric(format(first_week, "%Y"))) %>%
    count(first_week_year, first_week_num) %>%
    right_join(template_onsets) %>% ungroup()
    
onsets[is.na(onsets)] <- 0
onsets <- onsets %>%
    mutate(n = (n+1),
           sumn=sum(n),
           onset_prob = n/sum(n))

write.csv(onsets, file=paste0('inst/submissions/', filedate, '-onsets.csv'))



############################################
## calculuate next 4 week ahead bins      ##
############################################

ili_breaks <- seq(.5, 13, by = .5)
pred_bins <- preds_df %>%
    filter(week_date <= as.Date(paste0(last_obs_year, sprintf("%02d", last_obs_week), 00), 
                                format="%Y%W%w") + weeks(4),
           week_date > as.Date(paste0(last_obs_year, sprintf("%02d", last_obs_week), 00), 
                               format="%Y%W%w")) %>%
    mutate(ili_bin=cut(ili, breaks=ili_breaks, right=FALSE)) %>%
    count(week_date, ili_bin) %>%
    spread(week_date, n)

## NEED TO MAKE THIS BETTER
pred_bins[is.na(pred_bins)] <- 1
pred_bins_dodge <- rbind(rep(1, 5), 
    pred_bins,
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5), 
    rep(1, 5))
for(i in 2:ncol(pred_bins_dodge)){
    pred_bins_dodge[,i] <- pred_bins_dodge[,i]/sum(pred_bins_dodge[,i])
}
## need to make adjustment for zero predictions

write.csv(pred_bins_dodge, file=paste0('inst/submissions/', filedate, '-pred_bins.csv'))

## for point predictions 
preds_df %>%
    filter(week_date <= as.Date(paste0(last_obs_year, sprintf("%02d", last_obs_week), 00), format="%Y%W%w") + weeks(4),
           week_date > as.Date(paste0(last_obs_year, sprintf("%02d", last_obs_week), 00), format="%Y%W%w")) %>%
    group_by(week_date) %>%
    summarize(median(ili))


######################
## for peak height  ##
######################

peak_ht <- preds_df %>%
    group_by(sim) %>%
    summarize(peak_hts = max(ili)) %>%
    mutate(peak_ht=cut(peak_hts, breaks=ili_breaks, right=FALSE)) %>%
    count(peak_ht) 
peak_height_dodge <- rbind(rep(1, 2), rep(1, 2), rep(1, 2), rep(1, 2), rep(1, 2),
                           peak_ht,
                           rep(1, 2), rep(1, 2), rep(1, 2), rep(1, 2), rep(1, 2), rep(1, 2),
                           rep(1, 2), rep(1, 2), rep(1, 2), rep(1, 2), rep(1, 2), rep(1, 2), 
                           rep(1, 2), rep(1, 2), rep(1, 2))
peak_height_dodge[,2] <- peak_height_dodge[,2]/sum(peak_height_dodge[,2])

write.csv(peak_height_dodge, file=paste0('inst/submissions/', filedate, '-peak-height.csv'))

preds_df %>%
    group_by(sim) %>%
    summarize(peak_hts = max(ili)) %>%
    ungroup() %>% summarize(median(peak_hts))
    
    
############################################
## plot predictions sanity check          ##
############################################

preds_sum <- preds_df %>%
    group_by(week_date) %>%
    summarize(median_ili = median(ili),
              p05 = quantile(ili, .05),
              p95 = quantile(ili, .95))


ggplot(preds_sum, aes(x=week_date)) + 
    geom_line(aes(y=median_ili)) +
    geom_ribbon(aes(ymin=p05, ymax=p95), alpha=.2)

