
# data prep for east coast RCP analysis -----------------------------------

library(dplyr)
library(vegan)


floristics <- read.cep(file = "raw_floristic/PA_CANOCO.txt", maxdata = 240000000, positive = F, trace = T,  force = T)
floristics_site <- read.csv("raw_floristic/PA_ReplicateList.txt", stringsAsFactors = F, sep = "\t")
floristics$SiteNo <- floristics_site$SiteNo
save(floristics, file = "temp_PA_save.RData")
load("temp_PA_save.RData")

covariates <- read.csv("east_coast_covars.csv", stringsAsFactors = F)

# check out covariate data
cbind(names(covariates),
            unlist(lapply(X = 1:length(covariates), FUN = function(x, data){sum(is.na(data[,x]))}, covariates)))

# check out where the missing data is
plot(covariates$Latitude ~ covariates$Longitude, pch=16)
points(covariates$Latitude[is.na(covariates$cw_rain_sumwin_f)] ~ covariates$Longitude[is.na(covariates$cw_rain_sumwin_f)], col="red")
points(covariates$Latitude[is.na(covariates$finalsurplusAnn)] ~ covariates$Longitude[is.na(covariates$finalsurplusAnn)], col="red")
#plot(covariates$Latitude[covariates$SilicIndex!="No Value"] ~ covariates$Longitude[covariates$SilicIndex!="No Value"], pch=16, col=covariates$SilicIndex[covariates$SilicIndex!="No value"])
# ahh coastal effect!

#################################################
#################################################
# remmove missing data for the mean time
covariates <- covariates %>%
  filter(!is.na(cw_rain_sumwin_f), !is.na(finalsurplusAnn), !is.na(cw_prescott_f), SilicIndex!="No value")
# check
cbind(names(covariates),
      unlist(lapply(X = 1:length(covariates), FUN = function(x, data){sum(is.na(data[,x]))}, covariates)))
#################################################
#################################################


# COMBINE COVARS AND FLORISTICS -------------------------------------------

# remove what we don't need
covariates <- covariates %>%
  select(-(SilicName:SilicNA)) %>%
  mutate(SilicIndex = as.numeric(SilicIndex))
  
covariates.species <- inner_join(covariates, floristics, by="SiteNo")

# save off for model load
save(covariates.species, file = "covariates.species.RData")
