library(cdcfluview)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)

junk <- capture.output({
  usflu <- suppressMessages(get_flu_data("national", "ilinet", years=1997:2016))
})
ili_national <- suppressWarnings(transmute(usflu,
                                           region.type = REGION.TYPE,
                                           region = REGION,
                                           year = YEAR,
                                           week = WEEK,
                                           weighted_ili = as.numeric(X..WEIGHTED.ILI)))
ili_national$time <- ymd(paste(ili_national$year, "01", "01", sep = "-"))
week(ili_national$time) <- ili_national$week
ili_national$time_index <- seq_len(nrow(ili_national))

## Season column: for example, weeks of 2010 up through and including week 30 get season 2009/2010;
## weeks after week 30 get season 2010/2011
ili_national$season <- ifelse(
  ili_national$week <= 30,
  paste0(ili_national$year - 1, "/", ili_national$year),
  paste0(ili_national$year, "/", ili_national$year + 1)
)

## Season week column: week number within season
ili_national$season_week <- sapply(seq_len(nrow(ili_national)), function(row_ind) {
  sum(ili_national$season == ili_national$season[row_ind] &
        ili_national$time_index <= ili_national$time_index[row_ind])
})


## Subset to data actually used in this analysis -- up through end of 2014.
#ili_national <- ili_national[ili_national$year <= 2014, , drop = FALSE]

## cutoff time for training data
ili_train_cutoff_time <- ili_national$time[max(which(ili_national$year == 2010))]


christmas_week_inds <- sapply(
  unique(year(as.Date(ili_national$time))),
  function(year_val) {
    max(which(year(ili_national$time) == year_val &
                month(ili_national$time) == 12 &
                day(ili_national$time) <= 25))
  }
)

last_week_inds <- sapply(
  unique(year(as.Date(ili_national$time))),
  function(year_val) {
    max(which(year(ili_national$time) == year_val))
  }
)

## get legend grob
temp_df <- ili_national
temp_df$label <- "Incidence on\nChristmas Week"

p <- ggplot() +
  geom_point(aes(x = as.Date(time), y = weighted_ili, colour = label),
             size = 3,
             data = temp_df) +
  scale_colour_manual("", values = "grey") +
  theme_bw()

suppressWarnings(print(p))
grid_list <- capture.output(grid.ls())
legend_grob <- grid.get(str_trim(grid_list[grep("guide-box", grid_list)]))

panel_cutoff_ind <- 450
panel_cutoff_inds <- c(0, 300, 600, nrow(ili_national))


## set up layout for plot
pdf("/media/evan/data/Reich/ssr-influenza-competition/inst/wrapup-meetings/2015-2016/presentation/Plots/christmas-effect.pdf", width=14, height=6)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 3,
                                           heights = unit(c(1, 1, 1, 1), c("lines", "null", "null", "null")),
                                           widths = unit(c(1, 1, 0.3), c("lines", "null", "null")))))

## draw legend
pushViewport(viewport(layout.pos.row = 2:4, layout.pos.col = 3))
grid.draw(legend_grob)
upViewport()

#downViewport(xlab_grob_name)
grid.text("\`\`Christmas effect\'\' in Influenza-like Illness Incidence",
          #    x = unit(0.31, "npc"),
          #    y = unit(1, "npc"),
          #    rot = 90,
          gp = gpar(fontsize = 14),
          vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

grid.text("Weighted Influenza-like Illness",
          #    x = unit(0.31, "npc"),
          #    y = unit(1, "npc"),
          rot = 90,
          gp = gpar(fontsize = 12),
          vp = viewport(layout.pos.row = 2:4, layout.pos.col = 1))
#upViewport(0)


for(ind in 1:3) {
  p_ind <- ggplot() +
    geom_point(aes(x = as.Date(time), y = weighted_ili),
               size = 3,
               colour = "red",
               data = ili_national[christmas_week_inds[
                 christmas_week_inds > panel_cutoff_inds[ind] &
                   christmas_week_inds <= panel_cutoff_inds[ind + 1]], ]) +
    #        data = ili_national[last_week_inds, ]) +
    geom_line(aes(x = as.Date(time), y = weighted_ili),
              data = ili_national[seq(from = panel_cutoff_inds[ind] + 1, to = panel_cutoff_inds[ind + 1]), ]) +
    geom_vline(aes(xintercept = as.numeric(as.Date(time))),
               colour = "grey",
               data = ili_national[is.na(ili_national$weighted_ili) &
                                     seq_len(nrow(ili_national)) > panel_cutoff_inds[ind] &
                                     seq_len(nrow(ili_national)) <= panel_cutoff_inds[ind + 1], ]) +
    #    geom_vline(aes(xintercept = as.numeric(as.Date(ili_train_cutoff_time))),
    #        colour = "red", linetype = 2) +
    scale_x_date() +
    #    scale_x_date(limits = time_limits, expand = c(0, 0)) +
    xlab("Time") +
    ylab("") +
    ylim(c(0, 8)) +
    #    ylab("Weighted Influenza-like Illness\n") +
    #    ggtitle("Influenza Data - National United States") +
    theme_bw(base_size = 11)
  
  print(p_ind, vp = viewport(layout.pos.row = ind + 1, layout.pos.col = 2))
}

dev.off()
