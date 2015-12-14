## retrieve and save data for estimation
## Nicholas Reich
## December 2015

library(cdcfluview)
library(dplyr)

usflu<-get_flu_data("national", "ilinet", years=1997:2015)
usflu_to_save <- transmute(usflu,
                           region.type = REGION.TYPE,
                           region = REGION,
                           year = YEAR,
                           week = WEEK,
                           total_cases = X..WEIGHTED.ILI)
saveRDS(usflu_to_save, file='data/ili-national-1997-to-2015.rds')

for (i in 1:10){
    regionflu<-get_flu_data("HHS", sub_region=i, "ilinet", years=1997:2015)
    regionflu_to_save <- transmute(regionflu,
                                   region.type = REGION.TYPE,
                                   region = REGION,
                                   year = YEAR,
                                   week = WEEK,
                                   total_cases = X..WEIGHTED.ILI)
    filename <- paste0('data/ili-region', i, '-1997-to-2015.rds')
    saveRDS(regionflu_to_save, file=filename)
}