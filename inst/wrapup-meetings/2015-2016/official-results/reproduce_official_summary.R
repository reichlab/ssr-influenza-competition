library("dplyr")
library("tidyr")

setwd("/media/evan/data/Reich/ssr-influenza-competition/inst/wrapup-meetings/2015-2016/official-results")

results <- read.csv("KoT_Scores.csv", stringsAsFactors = FALSE)

relevant_results <- results %>%
  filter(location == "us", target == "onset")

sapply(seq_len(nrow(relevant_results)),
  function(i) {
    mean(relevant_results$score[seq_len(i)])
  })
