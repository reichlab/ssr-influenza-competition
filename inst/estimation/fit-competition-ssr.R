### get ssr fit -- prediction horizons 1:52
library(lubridate)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)
library(ssrFlu)
library(doMC)

registerDoMC(cores=3)

args <- commandArgs(trailingOnly=TRUE)

data_set <- args[1]
prediction_horizon_limit <- as.integer(args[2])

## load dataset
if(identical(data_set, "ili_national")) {
    data <- ili_national
} 
if(identical(data_set, "ili_region1")) {
    data <- ili_region1
} 
if(identical(data_set, "ili_region2")) {
    data <- ili_region2
} 
if(identical(data_set, "ili_region3")) {
    data <- ili_region3
} 
if(identical(data_set, "ili_region4")) {
    data <- ili_region4
} 
if(identical(data_set, "ili_region5")) {
    data <- ili_region5
} 
if(identical(data_set, "ili_region6")) {
    data <- ili_region6
} 
if(identical(data_set, "ili_region7")) {
    data <- ili_region7
} 
if(identical(data_set, "ili_region8")) {
    data <- ili_region8
} 
if(identical(data_set, "ili_region9")) {
    data <- ili_region9
} 
if(identical(data_set, "ili_region10")) {
    data <- ili_region10
} 



## add log column
data$log_total_cases <- log(data$total_cases + 1)

## add week_start_date
char_dates <- paste(data$year, data$week, "0")
data$week_start_date <- as.Date(char_dates, format="%Y %W %w")

## add smooth log column
sm <- loess(log_total_cases ~ as.numeric(week_start_date), data=data, span=12 / nrow(data))
data$smooth_log_cases <- sm$fitted

## add time column
data$time_ind <- seq_len(nrow(data))




ssr_control <- create_ssr_control(X_names=c("smooth_log_cases", "time_ind"),
    y_names="total_cases",
    time_name=NULL,
    max_lag=list(smooth_log_cases=1,
        time_ind=0),
    prediction_horizons=as.integer(prediction_horizon_limit),
    kernel_fns=list(smooth_log_cases="squared_exp_kernel",
        time_ind="periodic_kernel"),
    theta_est=list(smooth_log_cases="bw",
        time_ind="bw"),
    theta_fixed=list(time_ind=list(period=pi / 52)),
    theta_transform_fns=list(
        squared_exp_kernel=list(
            bw=list(transform="log",
                detransform="exp")
        ),
        periodic_kernel=list(
            bw=list(transform="log",
                detransform="exp")
        )
    ),
    crossval_buffer=52,
    loss_fn_name="mae_from_kernel_weights_and_centers",
    loss_fn_args=list())


ssr_fit <- ssr(X_names=c("smooth_log_cases", "time_ind"),
    y_names="total_cases",
    time_name=NULL,
    data=data,
    ssr_control=ssr_control)

save(ssr_fit,
    file=paste0("/home/ngr67a/2015-cdc-flu-competition/fit-competition-ssr-ph", prediction_horizon_limit, "-", data_set, ".Rdata"))
