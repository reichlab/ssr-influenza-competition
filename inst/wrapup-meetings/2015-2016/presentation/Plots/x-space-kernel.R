library(ssr)
library(ggplot2)

San_Juan_train$log_total_cases <- log(San_Juan_train$total_cases + 1)
sm <- loess(log_total_cases ~ as.numeric(week_start_date), data=San_Juan_train, span=12 / nrow(San_Juan_train))
San_Juan_train$smooth_log_cases <- sm$fitted
San_Juan_train$lag_1_smooth_log_cases <- lag(San_Juan_train$smooth_log_cases)

San_Juan_train <- San_Juan_train[-1, ]


ind_center <- which(San_Juan_train$season == "1995/1996" & San_Juan_train$season_week == 7)
center <- San_Juan_train[ind_center, c("lag_1_smooth_log_cases", "smooth_log_cases")]
San_Juan_train$dist_from_center <- sapply(seq_len(nrow(San_Juan_train)), function(row_ind) {
	sqrt(sum((San_Juan_train[row_ind, c("lag_1_smooth_log_cases", "smooth_log_cases")] - center)^2))
})
San_Juan_train$obs_weight <- exp(-1 * San_Juan_train$dist_from_center / 0.5)
San_Juan_train$obs_weight <- San_Juan_train$obs_weight / sum(San_Juan_train$obs_weight)

p <- ggplot() +
    geom_point(aes(x=lag_1_smooth_log_cases, y = smooth_log_cases), size = 48, colour="black", alpha = 0.05, data = San_Juan_train[ind_center, c("lag_1_smooth_log_cases", "smooth_log_cases"), drop = FALSE])+
    geom_point(aes(x=lag_1_smooth_log_cases, y = smooth_log_cases), size = 40, colour="black", alpha = 0.1, data = San_Juan_train[ind_center, c("lag_1_smooth_log_cases", "smooth_log_cases"), drop = FALSE])+
    geom_point(aes(x=lag_1_smooth_log_cases, y = smooth_log_cases), size = 32, colour="black", alpha = 0.2, data = San_Juan_train[ind_center, c("lag_1_smooth_log_cases", "smooth_log_cases"), drop = FALSE])+
    geom_point(aes(x=lag_1_smooth_log_cases, y = smooth_log_cases), size = 24, colour="black", alpha = 0.4, data = San_Juan_train[ind_center, c("lag_1_smooth_log_cases", "smooth_log_cases"), drop = FALSE])+
    geom_point(aes(x=lag_1_smooth_log_cases, y = smooth_log_cases), size = 16, colour="black", alpha = 0.6, data = San_Juan_train[ind_center, c("lag_1_smooth_log_cases", "smooth_log_cases"), drop = FALSE])+
    geom_path(aes(x=lag_1_smooth_log_cases, y = smooth_log_cases),
        colour = "grey",
        data=San_Juan_train[San_Juan_train$season %in% c("1990/1991", "1991/1992", "1992/1993", "1993/1994", "1994/1995") |
            (San_Juan_train$season == "1995/1996" & San_Juan_train$season_week <= 7), ]) +
    geom_point(aes(x=lag_1_smooth_log_cases, y = smooth_log_cases, colour = obs_weight),
        size = 3,
        data=San_Juan_train[San_Juan_train$season %in% c("1990/1991", "1991/1992", "1992/1993", "1993/1994", "1994/1995") |
            (San_Juan_train$season == "1995/1996" & San_Juan_train$season_week <= 7), ]) +
    scale_colour_gradient2("Time Point\nWeight", low = "#56B4E9", high = "#E69F00", midpoint = 0.002) +
#    scale_alpha_continuous("Observation Weight") +
#    scale_alpha_continuous("Observation Weight",
#        range = c(1, 0)) +
    xlab("Lag 1 Smoothed Log Cases") +
    ylab("Smoothed Log Cases") +
    ggtitle("Calculation of Time Point Weights by\nComputing Similarity of Lagged Observations") +
    theme_bw(base_size=22)

pdf("/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/plots/x-space-kernel.pdf", width=10, height=10)
print(p)
dev.off()

