### post-process ssr fits to obtain predictions
library(lubridate)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)
library(ssrFlu)
library(cdcfluview)

locations <- c("ili_national", paste0("ili_region", 1:10))
pred_hzns <- 1:30

##for(location in locations) {
for(location in "ili_national")
    ## collect fits for the given location
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
        
        ## subset is different than for estimation: adding 1 to index to account for rows removed for lag?
        data <- data[(262+1):nrow(data),]
        ## add time column
        data$time_ind <- seq_len(nrow(data))
    } else {
        data <- Iquitos_test
    }
    
    ## add time column
    ## data$time_ind <- seq_len(nrow(data))
    
    ## call function that outputs predictions to spreadsheets given ssr fits
    n_sims <- 10
    outfile_path <- "inst/competition-predictions"
    
    make_competition_forecasts(
        ssr_fits_by_prediction_horizon=ssr_fits_by_prediction_horizon_limit,
        n_sims=n_sims,
        data=data,
        outfile_path=outfile_path,
        location=location,
        last_obs_seasons = "2015/2016",
        last_obs_week = 49)
#}
