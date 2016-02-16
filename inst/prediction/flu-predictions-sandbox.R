### post-process ssr fits to obtain predictions
library(lubridate)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)
library(ssrFlu)
library(cdcfluview)
library(tidyr)

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

filedate <- '20160212'
last_obs_week <- 5
last_obs_year <- 2016
nsim <- 1000
pred_horizons <- 1:30
#pred_horizons <- 1:10  THIS ONE WORKS.
preds <- matrix(nrow=nsim*length(pred_horizons), ncol=3)
colnames(preds) <- c("hzn", "sim", "ili")
preds[,"hzn"] <- rep(pred_horizons, each=nsim)
preds[,"sim"] <- rep(1:nsim, times=length(pred_horizons))

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
    preds[idx,"ili"] <- simulate_from_weighted_kde(nsim, tmp)
}

upto_current_data_idx <- which(data$season=="2015/2016" & data$season_week==40):nrow(data)
data_sims <- data_frame(week = rep(data[upto_current_data_idx, "season_week"], each=nsim), 
                        sim = rep(1:nsim, times=length(upto_current_data_idx)),
                        season_week = 201500 + 100*(week<30) + week,
                        ili = rep(data[upto_current_data_idx, "total_cases"], each=nsim))

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
    rep(1, 5), 
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
    
