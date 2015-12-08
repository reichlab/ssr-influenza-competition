### functions to make predictions for the Dengue Forecasting competition

make_competition_forecasts <- function(ssr_fits_by_prediction_horizon,
    n_sims,
    data,
    outfile_path,
    location) {
    
    ## make a data frame with seasons and weeks at which we make predictions
    
    ## values used by competition organizers -- I think week 0 means week 52
    ## of the previous season.
    last_obs_seasons <- c("2005/2006",
        "2006/2007",
        "2007/2008",
        "2008/2009")
    
    last_obs_weeks <- seq(from = 0, to = 48, by = 4)
    
    results <- expand.grid(last_obs_week=last_obs_weeks,
        last_obs_season=last_obs_seasons,
        stringsAsFactors=FALSE)
    
    ## values we use -- translate week 0 to week 52 of the previous season
    results$last_obs_season_ind1 <- results$last_obs_season
    results$last_obs_week_ind1 <- results$last_obs_week
    results$last_obs_season_ind1[results$last_obs_week == 0] <-
        sapply(results$last_obs_season_ind1[results$last_obs_week == 0],
            function(season) {
                init_season <- as.numeric(substr(season, 1, 4))
                return(paste0(init_season - 1, "/", init_season))
            })
    results$last_obs_week_ind1[results$last_obs_week == 0] <- 52
    
    ## set row names
    rownames(results) <- paste0(results$last_obs_season,
        "_wk",
        results$last_obs_week)
    
    ## create a separate results data frame for each quantity we are predicting
    ## and add extra columns representing the quantities being predicted
    peak_incidence_results <- results
    peak_incidence_results$point <- NA
    
    peak_incidence_dist_cutoffs <- c(15 * seq(from = 0, to = 10), Inf)
    for(ind in seq_len(length(peak_incidence_dist_cutoffs) - 1)) {
        if(peak_incidence_dist_cutoffs[ind + 1] < Inf) {
            result_name <- paste0("p(",
                peak_incidence_dist_cutoffs[ind],
                "<=peak_incidence<",
                peak_incidence_dist_cutoffs[ind + 1],
                ")")
        } else {
            result_name <- paste0("p(",
                peak_incidence_dist_cutoffs[ind],
                "<=peak_incidence)")
        }
        peak_incidence_results[[result_name]] <- NA
    }
    
    peak_week_results <- results
    peak_week_results$point <- NA
    
    peak_week_dist_values <- seq_len(52)
    for(ind in seq_along(peak_week_dist_values)) {
        result_name <- paste0("p(peak_week=",
            peak_week_dist_values[ind],
            ")")
        peak_week_results[[result_name]] <- NA
    }
    
    season_incidence_results <- results
    season_incidence_results$point <- NA
    
    season_incidence_dist_cutoffs <- c(seq(from = 0, to = 1000, by = 100), Inf)
    for(ind in seq_len(length(season_incidence_dist_cutoffs) - 1)) {
        if(season_incidence_dist_cutoffs[ind + 1] < Inf) {
            result_name <- paste0("p(",
                season_incidence_dist_cutoffs[ind],
                "<=season_incidence<",
                season_incidence_dist_cutoffs[ind + 1],
                ")")
        } else {
            result_name <- paste0("p(",
                season_incidence_dist_cutoffs[ind],
                "<=season_incidence)")
        }
        season_incidence_results[[result_name]] <- NA
    }
    
    
    ## columns of each data frame where results are stored
    peak_incidence_results_prediction_cols <- seq(from=ncol(results) + 1,
        to=ncol(peak_incidence_results))
    peak_week_results_prediction_cols <- seq(from=ncol(results) + 1,
        to=ncol(peak_week_results))
    season_incidence_results_prediction_cols <- seq(from=ncol(results) + 1,
        to=ncol(season_incidence_results))
    
    
    ## get predictions for each combination of last_obs_season and last_obs_week
    for(results_row in seq_len(nrow(results))) {
        last_obs_week_ind0 <- results[results_row, "last_obs_week"]
        last_obs_week <- results[results_row, "last_obs_week_ind1"]
        last_obs_season <- results[results_row, "last_obs_season_ind1"]
        
        ## assemble weekly fits for the current season --
        ## either the observed counts if we are at or past that week, or
        ## the weights and kernel centers if we are not yet at that week.
        weekly_fits <- lapply(seq_len(52), function(prediction_week) {
            data_ind <- which(data$season == last_obs_season &
                data$season_week == prediction_week)
            if(prediction_week <= last_obs_week_ind0) {
                ## get observed value for the given week
                return(list(observed_value=data[data_ind, "total_cases"]))
            } else {
                ## get kernel weights and centers for prediction in the given
                ## week
                prediction_horizon <- prediction_week - last_obs_week_ind0
                ssr_fit <- ssr_fits_by_prediction_horizon[[prediction_horizon]]
                max_lag <- max(unlist(ssr_fit$lags_hat))
                prediction_data_inds <- seq(from=data_ind - max_lag,
                    to=data_ind)
                training_data_inds_to_drop <- seq(from=data_ind - 52,
                    to=data_ind + 52)
                training_data_inds_to_drop <- training_data_inds_to_drop[
                    training_data_inds_to_drop >= 1 &
                        training_data_inds_to_drop <= nrow(data)
                ]
                
                return(ssr_predict(ssr_fit,
                    prediction_data=data[prediction_data_inds, , drop=FALSE],
                    leading_rows_to_drop=max(unlist(ssr_fit$lags_hat)),
                    additional_training_rows_to_drop=training_data_inds_to_drop,
                    prediction_horizon=prediction_horizon,
                    normalize_weights=TRUE))
            }
        })
        
        ## get predictions
        results_one_season_week <- make_competition_forecasts_one_season_week(
            weekly_fits=weekly_fits,
            n_sims=n_sims)
        
        ## store predictions in corresponding results data frames
        peak_incidence_results[results_row,
            peak_incidence_results_prediction_cols] <- c(
                results_one_season_week$peak_incidence_pt_est,
                results_one_season_week$peak_incidence_dist_est)
        peak_week_results[results_row,
            peak_week_results_prediction_cols] <- c(
                results_one_season_week$peak_week_pt_est,
                results_one_season_week$peak_week_dist_est)
        season_incidence_results[results_row,
            season_incidence_results_prediction_cols] <- c(
                results_one_season_week$season_incidence_pt_est,
                results_one_season_week$season_incidence_dist_est)
    }
    
    ## transpose so that results are in format required for export
    results <- t(results)
    peak_incidence_results <- t(peak_incidence_results)
    peak_week_results <- t(peak_week_results)
    season_incidence_results <- t(season_incidence_results)
    
    ## save results data frames
    write.csv(peak_incidence_results[-(1:4), ],
        file=file.path(outfile_path, paste0("peakinc_", location, ".csv")))
    write.csv(peak_week_results[-(1:4), ],
        file=file.path(outfile_path, paste0("peakweek_", location, ".csv")))
    write.csv(season_incidence_results[-(1:4), ],
        file=file.path(outfile_path, paste0("seasoninc_", location, ".csv")))
}

make_competition_forecasts_by_trajectory <- function(ssr_fit,
    n_sims,
    data,
    outfile_path,
    location) {
    
    ## make a data frame with seasons and weeks at which we make predictions
    
    ## values used by competition organizers -- I think week 0 means week 52
    ## of the previous season.
    last_obs_seasons <- c("2005/2006",
        "2006/2007",
        "2007/2008",
        "2008/2009")
    
    last_obs_weeks <- seq(from = 0, to = 48, by = 4)
    
    results <- expand.grid(last_obs_week=last_obs_weeks,
        last_obs_season=last_obs_seasons,
        stringsAsFactors=FALSE)
    
    ## values we use -- translate week 0 to week 52 of the previous season
    results$last_obs_season_ind1 <- results$last_obs_season
    results$last_obs_week_ind1 <- results$last_obs_week
    results$last_obs_season_ind1[results$last_obs_week == 0] <-
        sapply(results$last_obs_season_ind1[results$last_obs_week == 0],
            function(season) {
                init_season <- as.numeric(substr(season, 1, 4))
                return(paste0(init_season - 1, "/", init_season))
            })
    results$last_obs_week_ind1[results$last_obs_week == 0] <- 52
    
    ## set row names
    rownames(results) <- paste0(results$last_obs_season,
        "_wk",
        results$last_obs_week)
    
    ## create a separate results data frame for each quantity we are predicting
    ## and add extra columns representing the quantities being predicted
    peak_incidence_results <- results
    peak_incidence_results$point <- NA
    
    peak_incidence_dist_cutoffs <- c(15 * seq(from = 0, to = 10), Inf)
    for(ind in seq_len(length(peak_incidence_dist_cutoffs) - 1)) {
        if(peak_incidence_dist_cutoffs[ind + 1] < Inf) {
            result_name <- paste0("p(",
                peak_incidence_dist_cutoffs[ind],
                "<=peak_incidence<",
                peak_incidence_dist_cutoffs[ind + 1],
                ")")
        } else {
            result_name <- paste0("p(",
                peak_incidence_dist_cutoffs[ind],
                "<=peak_incidence)")
        }
        peak_incidence_results[[result_name]] <- NA
    }
    
    peak_week_results <- results
    peak_week_results$point <- NA
    
    peak_week_dist_values <- seq_len(52)
    for(ind in seq_along(peak_week_dist_values)) {
        result_name <- paste0("p(peak_week=",
            peak_week_dist_values[ind],
            ")")
        peak_week_results[[result_name]] <- NA
    }
    
    season_incidence_results <- results
    season_incidence_results$point <- NA
    
    season_incidence_dist_cutoffs <- c(seq(from = 0, to = 1000, by = 100), Inf)
    for(ind in seq_len(length(season_incidence_dist_cutoffs) - 1)) {
        if(season_incidence_dist_cutoffs[ind + 1] < Inf) {
            result_name <- paste0("p(",
                season_incidence_dist_cutoffs[ind],
                "<=season_incidence<",
                season_incidence_dist_cutoffs[ind + 1],
                ")")
        } else {
            result_name <- paste0("p(",
                season_incidence_dist_cutoffs[ind],
                "<=season_incidence)")
        }
        season_incidence_results[[result_name]] <- NA
    }
    
    
    ## columns of each data frame where results are stored
    peak_incidence_results_prediction_cols <- seq(from=ncol(results) + 1,
        to=ncol(peak_incidence_results))
    peak_week_results_prediction_cols <- seq(from=ncol(results) + 1,
        to=ncol(peak_week_results))
    season_incidence_results_prediction_cols <- seq(from=ncol(results) + 1,
        to=ncol(season_incidence_results))
    
    
    ## get predictions for each combination of last_obs_season and last_obs_week
    for(results_row in seq_len(nrow(results))) {
        last_obs_week_ind0 <- results[results_row, "last_obs_week"]
        last_obs_week <- results[results_row, "last_obs_week_ind1"]
        last_obs_season <- results[results_row, "last_obs_season_ind1"]
        
        if(last_obs_season == "2006/2007" & last_obs_week == "44") {
            browser()
        }
        
        ## get inds for prediction and training data -- ensure they're the
        ## same for all weeks in the season so that we can perform prediction
        ## by trajectory
        
        ## prediction data are those within max_lag of the last observed week
        max_lag <- max(unlist(ssr_fit$lags_hat))
        last_obs_week_ind <- which(data$season == last_obs_season &
                data$season_week == last_obs_week)
        prediction_data_inds <- seq(
            from=last_obs_week_ind - max_lag,
            to=last_obs_week_ind)
        
        ## drop indices that are within 1 year of the last observed week or
        ## within the last (prediction_horizon = 52 - last_obs_week_ind0)
        ## weeks of the end of the data
        training_data_inds_to_drop <- c(
            seq(from=last_obs_week_ind - 52,
                to=last_obs_week_ind + 52),
            seq(from=nrow(data) - (52 - last_obs_week_ind0),
                to=nrow(data)))
        training_data_inds_to_drop <- training_data_inds_to_drop[
            training_data_inds_to_drop >= 1 &
                training_data_inds_to_drop <= nrow(data)
            ]
        
        ## assemble weekly fits for the current season --
        ## either the observed counts if we are at or past that week, or
        ## the weights and kernel centers if we are not yet at that week.
        weekly_fits <- lapply(seq_len(52), function(prediction_week) {
            data_ind <- which(data$season == last_obs_season &
                    data$season_week == prediction_week)
            if(prediction_week <= last_obs_week_ind0) {
                ## get observed value for the given week
                return(list(observed_value=data[data_ind, "total_cases"]))
            } else {
                ## get kernel weights and centers for prediction in the given
                ## week
                prediction_horizon <- prediction_week - last_obs_week_ind0
                
                return(ssr_predict(ssr_fit,
                    prediction_data=data[prediction_data_inds, , drop=FALSE],
                    leading_rows_to_drop=max(unlist(ssr_fit$lags_hat)),
                    additional_training_rows_to_drop=training_data_inds_to_drop,
                    prediction_horizon=prediction_horizon,
                    normalize_weights=TRUE))
            }
        })
        
        ## get predictions
        results_one_season_week <-
            make_competition_forecasts_one_season_week_by_trajectory(
                weekly_fits=weekly_fits,
                n_sims=n_sims)
        
        ## store predictions in corresponding results data frames
        peak_incidence_results[results_row,
            peak_incidence_results_prediction_cols] <- c(
                results_one_season_week$peak_incidence_pt_est,
                results_one_season_week$peak_incidence_dist_est)
        peak_week_results[results_row,
            peak_week_results_prediction_cols] <- c(
                results_one_season_week$peak_week_pt_est,
                results_one_season_week$peak_week_dist_est)
        season_incidence_results[results_row,
            season_incidence_results_prediction_cols] <- c(
                results_one_season_week$season_incidence_pt_est,
                results_one_season_week$season_incidence_dist_est)
    }
    
    ## transpose so that results are in format required for export
    results <- t(results)
    peak_incidence_results <- t(peak_incidence_results)
    peak_week_results <- t(peak_week_results)
    season_incidence_results <- t(season_incidence_results)
    
    ## save results data frames
    write.csv(peak_incidence_results[-(1:4), ],
        file=file.path(outfile_path, paste0("peakinc_", location, ".csv")))
    write.csv(peak_week_results[-(1:4), ],
        file=file.path(outfile_path, paste0("peakweek_", location, ".csv")))
    write.csv(season_incidence_results[-(1:4), ],
        file=file.path(outfile_path, paste0("seasoninc_", location, ".csv")))
}

make_competition_forecasts_one_season_week <- function(weekly_fits,
    n_sims) {
    
    ## simulate n_sims realizations from the distributions implied by the
    ## weekly_kernel_fits
    ## result is a data frame with simulation index in the row,
    ## week of the season in the column
    simulated_counts_by_week <-
        simulate_counts_by_week(weekly_fits, n_sims)
    
    ## extract estimates of peak incidence
    max_counts_by_sim_ind <- apply(simulated_counts_by_week, 1, max)
    
    peak_incidence_pt_est <- mean(max_counts_by_sim_ind)
    
    peak_incidence_dist_cutoffs <- c(15 * seq(from = 0, to = 10), Inf)
    peak_incidence_dist_est <- sapply(
        seq_len(length(peak_incidence_dist_cutoffs) - 1),
        function(lb_ind) {
            mean(max_counts_by_sim_ind >= peak_incidence_dist_cutoffs[lb_ind] &
                max_counts_by_sim_ind < peak_incidence_dist_cutoffs[lb_ind + 1])
        })
    
    ## extract estimates of peak week
    peak_week_by_sim_ind <- apply(simulated_counts_by_week, 1, which.max)
    
    peak_week_pt_est <- mean(peak_week_by_sim_ind)
    
    peak_week_dist_est <- sapply(seq_len(52), function(week) {
        mean(peak_week_by_sim_ind == week)
    })
    
    ## extract estimates of season incidence
    season_incidence_by_sim_ind <- apply(simulated_counts_by_week, 1, sum)
    
    season_incidence_pt_est <- mean(season_incidence_by_sim_ind)
    
    season_incidence_dist_cutoffs <- c(seq(from = 0, to = 1000, by = 100), Inf)
    season_incidence_dist_est <- sapply(
        seq_len(length(season_incidence_dist_cutoffs) - 1),
        function(lb_ind) {
            mean(season_incidence_by_sim_ind >= season_incidence_dist_cutoffs[lb_ind] &
                season_incidence_by_sim_ind < season_incidence_dist_cutoffs[lb_ind + 1])
        })
    
    return(list(
        peak_incidence_pt_est=peak_incidence_pt_est,
        peak_incidence_dist_est=peak_incidence_dist_est,
        peak_week_pt_est=peak_week_pt_est,
        peak_week_dist_est=peak_week_dist_est,
        season_incidence_pt_est=season_incidence_pt_est,
        season_incidence_dist_est=season_incidence_dist_est
    ))
}



make_competition_forecasts_one_season_week_by_trajectory <-
    function(weekly_fits,
        n_sims) {
    
    ## simulate n_sims realizations from the distributions implied by the
    ## weekly_kernel_fits
    ## result is a data frame with simulation index in the row,
    ## week of the season in the column
    simulated_counts_by_week <-
        simulate_counts_by_week_by_trajectory(weekly_fits, n_sims)
    
    ## extract estimates of peak incidence
    max_counts_by_sim_ind <- apply(simulated_counts_by_week, 1, max)
    
    peak_incidence_pt_est <- mean(max_counts_by_sim_ind)
    
    peak_incidence_dist_cutoffs <- c(15 * seq(from = 0, to = 10), Inf)
    peak_incidence_dist_est <- sapply(
        seq_len(length(peak_incidence_dist_cutoffs) - 1),
        function(lb_ind) {
            mean(max_counts_by_sim_ind >= peak_incidence_dist_cutoffs[lb_ind] &
                    max_counts_by_sim_ind < peak_incidence_dist_cutoffs[lb_ind + 1])
        })
    
    ## extract estimates of peak week
    peak_week_by_sim_ind <- apply(simulated_counts_by_week, 1, which.max)
    
    peak_week_pt_est <- mean(peak_week_by_sim_ind)
    
    peak_week_dist_est <- sapply(seq_len(52), function(week) {
        mean(peak_week_by_sim_ind == week)
    })
    
    ## extract estimates of season incidence
    season_incidence_by_sim_ind <- apply(simulated_counts_by_week, 1, sum)
    
    season_incidence_pt_est <- mean(season_incidence_by_sim_ind)
    
    season_incidence_dist_cutoffs <- c(seq(from = 0, to = 1000, by = 100), Inf)
    season_incidence_dist_est <- sapply(
        seq_len(length(season_incidence_dist_cutoffs) - 1),
        function(lb_ind) {
            mean(season_incidence_by_sim_ind >= season_incidence_dist_cutoffs[lb_ind] &
                    season_incidence_by_sim_ind < season_incidence_dist_cutoffs[lb_ind + 1])
        })
    
    return(list(
        peak_incidence_pt_est=peak_incidence_pt_est,
        peak_incidence_dist_est=peak_incidence_dist_est,
        peak_week_pt_est=peak_week_pt_est,
        peak_week_dist_est=peak_week_dist_est,
        season_incidence_pt_est=season_incidence_pt_est,
        season_incidence_dist_est=season_incidence_dist_est
    ))
}

simulate_counts_by_week <- function(weekly_fits, n_sims) {
    results <- sapply(seq_len(52), function(week) {
        if(length(weekly_fits[[week]]$observed_value) == 1) {
            return(rep(weekly_fits[[week]]$observed_value, n_sims))
        } else {
            return(simulate_from_weighted_kde(n_sims, weekly_fits[[week]]))
        }
    })
    
    return(results)
}

simulate_counts_by_week_by_trajectory <- function(weekly_fits, n_sims) {
    ## get weights -- same for all weekly fits.
    weights <- weekly_fits[[52]]$weights
    
    ## select index of the trajectory to sample from for each simulation
    trajectory_index <- sample(length(weights),
        size=n_sims,
        replace=TRUE,
        prob=weights)
    
    ## return observed value if available or simulated value
    results <- sapply(seq_len(52), function(week) {
        if(length(weekly_fits[[week]]$observed_value) == 1) {
            return(rep(weekly_fits[[week]]$observed_value, n_sims))
        } else {
            return(simulate_from_weighted_kde_given_ind(trajectory_index,
                weekly_fits[[week]]))
        }
    })

    return(results)
}

simulate_from_weighted_kde_given_ind <- function(inds, weighted_kde_fit) {
    ## get bandwidth
    density_fit <- density(x=weighted_kde_fit$centers[, 1],
        weights=weighted_kde_fit$weights,
        bw="SJ")
    
    ## get simulated values -- from a mixture of normals with sd=density_fit$bw
    component_means <- weighted_kde_fit$centers[inds, 1]
    
    return(rnorm(length(inds), component_means, density_fit$bw))
}

simulate_from_weighted_kde <- function(n, weighted_kde_fit) {
    ## get bandwidth
    density_fit <- density(x=weighted_kde_fit$centers[, 1],
        weights=weighted_kde_fit$weights,
        bw="SJ")
    
    ## get simulated values -- from a mixture of normals with sd=density_fit$bw
    component_means <- sample(weighted_kde_fit$centers[, 1],
        size=n,
        replace=TRUE,
        prob=weighted_kde_fit$weights)
    
    return(rnorm(n, component_means, density_fit$bw))
}
