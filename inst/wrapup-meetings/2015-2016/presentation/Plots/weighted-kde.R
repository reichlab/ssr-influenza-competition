library(ggplot2)

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

p <- ggplot() +
	geom_point(aes(x = x, y = y, colour = col),
		data = data.frame(x = x, y = rep(0, 5), col = factor("dummy_col", levels = c("dummy_col", "blah"))),
		size = 5) +
	geom_line(aes(x = x, y = y), linetype = 2, data = kde_df_1) +
	geom_line(aes(x = x, y = y), linetype = 2, data = kde_df_2) +
	geom_line(aes(x = x, y = y), linetype = 2, data = kde_df_3) +
	geom_line(aes(x = x, y = y), linetype = 2, data = kde_df_4) +
	geom_line(aes(x = x, y = y), linetype = 2, data = kde_df_5) +
	geom_line(aes(x = x, y = y), linetype = 1, data = kde_df_combined) +
	scale_colour_manual("Points",
		breaks = c("Observations", ""),
		values = "#E69F00") +
	xlim(c(0, 100)) +
	ylim(c(0, 0.05)) +
	xlab("Y") +
	ylab("Estimated Density") +
	ggtitle("Weighted Kernel Density Estimation") +
	theme_bw(base_size = 22) +
	theme()


p <- ggplot() +
	geom_line(aes(x = x, y = y, linetype = curve_type, group = curve_group), colour = "blue", data = plot_df) +
	geom_point(aes(x = x, y = y, colour = col),
		data = data.frame(x = x, y = rep(0, 5), col = factor("dummy_col", levels = c("dummy_col", "blah"))),
		size = 5) +
	scale_colour_manual("Points",
		labels = "Observed Value\nof Y",
		values = "#E69F00") +
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


pdf("/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/plots/conditional-kde-overview.pdf", width=10, height=4)
print(p)
dev.off()

