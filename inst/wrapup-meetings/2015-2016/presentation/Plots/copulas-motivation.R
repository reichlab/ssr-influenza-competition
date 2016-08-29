library(ggplot2)

grid_length <- 100L
grid_vals <- seq(from = 0, to = 10, length = grid_length)
plot_df <- data.frame(
  horizon = factor(c(rep(1, grid_length), rep(2, grid_length), rep(3, grid_length), rep(4, grid_length), rep(5, grid_length))),
  y = c(grid_vals, grid_vals, grid_vals, grid_vals, grid_vals),
  x = c(1 - dnorm(grid_vals, mean = 5, sd = 0.8),
        2 - dnorm(grid_vals, mean = 5, sd = 0.8),
        3 - dnorm(grid_vals, mean = 5, sd = 0.8),
        4 - dnorm(grid_vals, mean = 5, sd = 0.8),
        5 - dnorm(grid_vals, mean = 5, sd = 0.8))
)

possible_traj_df <- data.frame(
  x = c(1:5, 1:5),
  y = c(rep(6, 5), 6, 4, 6, 4, 6),
  traj_ind = factor(c(rep(1, 5), rep(2, 5)))
)

p <- ggplot() +
  geom_path(aes(x = x, y = y, group = horizon), data = plot_df) +
  geom_path(aes(x = x, y = y, group = traj_ind, colour = traj_ind, linetype = traj_ind), size = 1, data = possible_traj_df) +
  xlab("Prediction Horizon") +
  ylab("Incidence") +
  scale_x_continuous(breaks = 1:5) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")

pdf("/media/evan/data/Reich/ssr-influenza-competition/inst/wrapup-meetings/2015-2016/presentation/Plots/copulas-motivation.pdf", width=10, height=4)
print(p)
dev.off()
