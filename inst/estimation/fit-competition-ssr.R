### get ssr fit -- prediction horizons 1:52
library(lubridate)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)
library(ssr)
library(doMC)

registerDoMC(cores=3)

args <- commandArgs(trailingOnly=TRUE)

data_set <- args[1]
prediction_horizon_limit <- as.integer(args[2])

if(identical(data_set, "San_Juan")) {
    data <- San_Juan_train
} else {
    data <- Iquitos_train
}

## add log column
data$log_total_cases <- log(data$total_cases + 1)

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
    file=paste0("/home/er71a/ssr-poster/fit-competition-ssr-ph", prediction_horizon_limit, "-", data_set, ".Rdata"))
