### post-process ssr fits to obtain predictions
library(lubridate)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)
library(ssr)

for(location in c("sanjuan", "iquitos")) {
#for(location in c("iquitos")) {
#for(location in c("sanjuan")) {
    ## collect fits for the given location
    if(identical(location, "sanjuan")) {
        location_for_ssr_fit_file <- "San_Juan"
    } else {
        location_for_ssr_fit_file <- "Iquitos"
    }
    
    ssr_fits_by_prediction_horizon_limit <- lapply(seq(from=4, to=52, by=4),
        function(phl) {
            file_name <- paste0(
                "F:/Reich/dengue-ssr-prediction/competition-fits/fit-competition-ssr-ph",
                phl,
                "-",
                location_for_ssr_fit_file,
                ".Rdata")
            read_env <- new.env()
            load(file_name, envir=read_env)
            return(read_env$ssr_fit)
        })
    names(ssr_fits_by_prediction_horizon_limit) <-
        paste0("phl", seq(from=4, to=52, by=4))
    
    ## assemble data set
    if(identical(location, "sanjuan")) {
        data <- San_Juan_test
    } else {
        data <- Iquitos_test
    }
    
    ## add time column
    data$time_ind <- seq_len(nrow(data))
    
    ## call function that outputs predictions to spreadsheets given ssr fits
    n_sims <- 100000
    outfile_path <- "F:/Reich/dengue-ssr-prediction/competition-predictions"
    
    make_competition_forecasts_by_trajectory(
        ssr_fits_by_prediction_horizon_limit=ssr_fits_by_prediction_horizon_limit,
        n_sims=n_sims,
        data=data,
        outfile_path=outfile_path,
        location=location,
        phase="test")
}
