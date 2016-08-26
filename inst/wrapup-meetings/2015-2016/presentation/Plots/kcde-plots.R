suppressWarnings(library(knitr))
suppressWarnings(library(plyr))
suppressWarnings(library(grid))
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))

San_Juan_test <- read.csv("/media/evan/data/Reich/ssr-influenza-competition/inst/wrapup-meetings/2015-2016/presentation/Plots/San_Juan_Testing_Data.csv")
San_Juan_test$week_start_date <- as.Date(San_Juan_test$week_start_date)

subset_inds <- which(San_Juan_test$season %in% c("2005/2006", "2006/2007", "2007/2008", "2008/2009", "2009/2010", "2010/2011", "2011/2012", "2012/2013") &
                       !(San_Juan_test$season == "2012/2013" & San_Juan_test$season_week > 18))
San_Juan_subset <- San_Juan_test[subset_inds, ]

recent_inds <- seq(from=nrow(San_Juan_subset) - 1,
                   to = nrow(San_Juan_subset))

similar_inds <- c(14, 15, 121, 122, 265, 266)









subset_inds <- which(San_Juan_test$season %in% c("2005/2006", "2006/2007", "2007/2008", "2008/2009", "2009/2010", "2010/2011", "2011/2012", "2012/2013") &
                       !(San_Juan_test$season == "2012/2013" & San_Juan_test$season_week > 18))
San_Juan_subset <- San_Juan_test[subset_inds, ]

recent_inds <- seq(from=nrow(San_Juan_subset) - 1,
                   to = nrow(San_Juan_subset))
#similar_inds <- which(abs(San_Juan_subset$total_cases - San_Juan_subset$total_cases[nrow(San_Juan_subset)]) < 15)
similar_inds <- c(14, 15, 121, 122, 265, 266)
prediction_length <- 4
prediction_inds <- c(15 + prediction_length,
                     122 + prediction_length,
                     266 + prediction_length)

error_sds <- c(3, 4, 7, 9, 13, 15, 12, 11, 10, 9, 9, 8)
sd_multipliers <- c(0.1, 1, 2, 4)

dist_means_df <- data.frame(
  ph = prediction_length,
  est_mean = mean(San_Juan_subset[c(15, 122, 266) + prediction_length, "total_cases"])
)

alpha_levels <- c(1, 0.8, 0.6, 0.4)

dist_df <- rbind.fill(lapply(prediction_length, function(ph) {
  est_mean <- dist_means_df$est_mean[ph]
  data.frame(
    week_start_date=San_Juan_subset[nrow(San_Juan_subset), "week_start_date"] + ph,
    ph = ph,
    alpha_level = alpha_levels,
    y_min = est_mean - sd_multipliers * error_sds[ph],
    y_max = est_mean + sd_multipliers * error_sds[ph]
  )
}))

pt_size <- 4
p <- ggplot(San_Juan_subset) +
  geom_line(aes(x = week_start_date, y = total_cases)) +
  # geom_linerange(aes(x = week_start_date, ymin = y_min, ymax = y_max, alpha = alpha_level),
  #     colour = "blue",
  #     size = pt_size,
  #     data=dist_df) +
  geom_point(aes(x = week_start_date, y = total_cases),
             colour = "blue",
             size = pt_size,
             shape = 16,
             data=data.frame(week_start_date=San_Juan_subset[similar_inds, "week_start_date"] + prediction_length,
                             total_cases=San_Juan_subset[similar_inds, "total_cases"])) +
  geom_point(aes(x = week_start_date, y = total_cases),
             colour = "blue",
             size = pt_size,
             shape = 15,
             data=data.frame(week_start_date=San_Juan_subset[c(15, 122, 266), "week_start_date"] + prediction_length * 7,
                             total_cases=San_Juan_subset[c(15, 122, 266) + prediction_length, "total_cases"])) +
  geom_point(aes(x = week_start_date, y = total_cases),
             colour = "orange",
             size = pt_size,
             shape = 18,
             data=data.frame(week_start_date=San_Juan_subset[c(nrow(San_Juan_subset) - 1, nrow(San_Juan_subset)), "week_start_date"] + prediction_length,
                             total_cases=San_Juan_subset[c(nrow(San_Juan_subset) - 1, nrow(San_Juan_subset)), "total_cases"])) +
  geom_point(aes(x = week_start_date, y = total_cases),
             colour = "orange",
             size = pt_size,
             shape = 17,
             data=data.frame(week_start_date=San_Juan_subset[nrow(San_Juan_subset), "week_start_date"] + prediction_length * 7,
                             total_cases=sapply(prediction_length, function(ph) {
                               mean(San_Juan_subset[c(15, 122, 266) + ph, "total_cases"])
                             }))) +
  xlab("Time") +
  ylab("Weekly\nCases") +
  theme_bw(base_size = 20)

png("/media/evan/data/Reich/ssr-influenza-competition/inst/wrapup-meetings/2015-2016/presentation/Plots/moa-overview.png", width=850, height=250)
print(p)
dev.off()


















San_Juan_train <- read.csv("/media/evan/data/Reich/ssr-influenza-competition/inst/wrapup-meetings/2015-2016/presentation/Plots/San_Juan_Training_Data.csv")
San_Juan_train$week_start_date <- as.Date(San_Juan_train$week_start_date)

San_Juan_train$log_total_cases <- log(San_Juan_train$total_cases + 1)
sm <- loess(log_total_cases ~ as.numeric(week_start_date), data=San_Juan_train, span=12 / nrow(San_Juan_train))
San_Juan_train$smooth_log_cases <- sm$fitted
San_Juan_train$lag_1_smooth_log_cases <- c(NA, San_Juan_train$smooth_log_cases[seq_len(nrow(San_Juan_train) - 1)])

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
             size = 5,
             data=San_Juan_train[San_Juan_train$season %in% c("1990/1991", "1991/1992", "1992/1993", "1993/1994", "1994/1995") |
                                   (San_Juan_train$season == "1995/1996" & San_Juan_train$season_week <= 7), ]) +
  geom_point(aes(x = lag_1_smooth_log_cases, y = smooth_log_cases),
             colour = "orange",
             shape = 18,
             size = 6,
             data = San_Juan_train[ind_center, c("lag_1_smooth_log_cases", "smooth_log_cases"), drop = FALSE]) +
  scale_colour_gradient2("Time\nPoint\nWeight", low = "gray", mid = "lightblue", high = "blue", midpoint = 0.002) +
  #    scale_alpha_continuous("Observation Weight") +
  #    scale_alpha_continuous("Observation Weight",
  #        range = c(1, 0)) +
  xlab("Lag 1 Smoothed Log Cases") +
  ylab("Smoothed Log Cases") +
  ggtitle("Calculation of Time Point Weights by\nComputing Similarity of Lagged Observations") +
  theme_bw(base_size=22)

pdf("/media/evan/data/Reich/ssr-influenza-competition/inst/wrapup-meetings/2015-2016/presentation/Plots/x-space-kernel.pdf", width=10, height=10)
print(p)
dev.off()































x <- c(17, 23, 56, 65, 71)
bw <- 5

x_grid <- seq(from = 0, to = 100)



kde_df_1 <- data.frame(
  x = x_grid,
  y = dnorm(x_grid, mean = x[1], sd = bw) * 0.05,
  curve_group = "k1",
  curve_type = "weighted_density"
)
kde_df_2 <- data.frame(
  x = x_grid,
  y = dnorm(x_grid, mean = x[2], sd = bw) * 0.07,
  curve_group = "k2",
  curve_type = "weighted_density"
)
kde_df_3 <- data.frame(
  x = x_grid,
  y = dnorm(x_grid, mean = x[3], sd = bw) * 0.2,
  curve_group = "k3",
  curve_type = "weighted_density"
)
kde_df_4 <- data.frame(
  x = x_grid,
  y = dnorm(x_grid, mean = x[4], sd = bw) * 0.4,
  curve_group = "k4",
  curve_type = "weighted_density"
)
kde_df_5 <- data.frame(
  x = x_grid,
  y = dnorm(x_grid, mean = x[5], sd = bw) * 0.28,
  curve_group = "k5",
  curve_type = "weighted_density"
)

kde_df_combined <- data.frame(
  x = x_grid,
  y = apply(cbind(kde_df_1$y, kde_df_2$y, kde_df_3$y, kde_df_4$y, kde_df_5$y),
            1,
            sum),
  curve_group = "kde",
  curve_type = "estimated_density"
)

plot_df <- rbind.fill(kde_df_1, kde_df_2, kde_df_3, kde_df_4, kde_df_5, kde_df_combined)

# p <- ggplot() +
#   geom_point(aes(x = x, y = y, colour = col),
#              data = data.frame(x = x, y = rep(0, 5), col = factor("dummy_col", levels = c("dummy_col", "blah"))),
#              size = 5) +
#   geom_line(aes(x = x, y = y), linetype = 2, data = kde_df_1) +
#   geom_line(aes(x = x, y = y), linetype = 2, data = kde_df_2) +
#   geom_line(aes(x = x, y = y), linetype = 2, data = kde_df_3) +
#   geom_line(aes(x = x, y = y), linetype = 2, data = kde_df_4) +
#   geom_line(aes(x = x, y = y), linetype = 2, data = kde_df_5) +
#   geom_line(aes(x = x, y = y), linetype = 1, data = kde_df_combined) +
#   scale_colour_manual("Points",
#                       breaks = c("Observations", ""),
#                       values = "orange", "") +
#   xlim(c(0, 100)) +
#   ylim(c(0, 0.05)) +
#   xlab("Y") +
#   ylab("Estimated Density") +
#   ggtitle("Weighted Kernel Density Estimation") +
#   theme_bw(base_size = 22) +
#   theme()


p <- ggplot() +
  geom_line(aes(x = x, y = y, linetype = curve_type, group = curve_group), colour = "orange", data = plot_df) +
  geom_point(aes(x = x, y = y, colour = col),
             shape = 15,
             data = data.frame(x = x, y = rep(0, 5), col = factor("dummy_col", levels = c("dummy_col", "blah"))),
             size = 5) +
  scale_colour_manual("Points",
                      labels = "Observed Value\nof Y",
                      values = "blue") +
  scale_linetype_manual("Curves",
                        breaks = c("weighted_density", "estimated_density"),
                        labels = c("Weighted\nKernel Function", "Combined\nDensity Estimate"),
                        values = c(2, 1)) +
  xlim(c(0, 100)) +
  ylim(c(0, 0.05)) +
  xlab("Y") +
  ylab("Estimated Density") +
  ggtitle("Weighted Kernel Density Estimation") +
  theme_bw(base_size = 22) +
  theme()


pdf("/media/evan/data/Reich/ssr-influenza-competition/inst/wrapup-meetings/2015-2016/presentation/Plots/weighted-kde.pdf", width=10, height=7)
print(p)
dev.off()
