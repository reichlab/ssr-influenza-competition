## retrieve and save data for estimation
## Nicholas Reich
## December 2015

library(cdcfluview)
library(dplyr)

usflu<-get_flu_data("national", "ilinet", years=1997:2015)
ili_national <- transmute(usflu,
                           region.type = REGION.TYPE,
                           region = REGION,
                           year = YEAR,
                           week = WEEK,
                           total_cases = as.numeric(X..WEIGHTED.ILI))
save(ili_national, file='data/ili-national-1997-to-2015.rda')

regionflu<-get_flu_data("HHS", sub_region=1:10, "ilinet", years=1997:2015)
regionflu_to_save <- transmute(regionflu,
                               region.type = REGION.TYPE,
                               region = REGION,
                               year = YEAR,
                               week = WEEK,
                               total_cases = as.numeric(X..WEIGHTED.ILI))

## save each region separately
ili_region1 <- dplyr::filter(regionflu_to_save, region == "Region 1")
save(ili_region1, file='data/ili-region1-1997-to-2015.rda')

ili_region2 <- dplyr::filter(regionflu_to_save, region == "Region 2")
save(ili_region2, file='data/ili-region2-1997-to-2015.rda')

ili_region3 <- dplyr::filter(regionflu_to_save, region == "Region 3")
save(ili_region3, file='data/ili-region3-1997-to-2015.rda')

ili_region4 <- dplyr::filter(regionflu_to_save, region == "Region 4")
save(ili_region4, file='data/ili-region4-1997-to-2015.rda')

ili_region5 <- dplyr::filter(regionflu_to_save, region == "Region 5")
save(ili_region5, file='data/ili-region5-1997-to-2015.rda')

ili_region6 <- dplyr::filter(regionflu_to_save, region == "Region 6")
save(ili_region6, file='data/ili-region6-1997-to-2015.rda')

ili_region7 <- dplyr::filter(regionflu_to_save, region == "Region 7")
save(ili_region7, file='data/ili-region7-1997-to-2015.rda')

ili_region8 <- dplyr::filter(regionflu_to_save, region == "Region 8")
save(ili_region8, file='data/ili-region8-1997-to-2015.rda')

ili_region9 <- dplyr::filter(regionflu_to_save, region == "Region 9")
save(ili_region9, file='data/ili-region9-1997-to-2015.rda')

ili_region10 <- dplyr::filter(regionflu_to_save, region == "Region 10")
save(ili_region10, file='data/ili-region10-1997-to-2015.rda')

