library(methods) # workaround so functions can be found when running via RScript
library(Matrix) # sameo
library(glmnet) # sameo
library(RCPmod) # normally, this is the only package you'll need to load, rest will load via namespace

# data prep ---------------------------------------------------------------
## load data that is prepped from "east_coast_data.R"
## alternatively, load or create a single object that is site by species/covariates

# covariates.species = read.csv("covariates_species.csv", 
#                               header=TRUE, stringsAsFactors=FALSE)

# loading an .RData file, since uncompressed .csv is a little unwieldly in this case
load("covariates.species.RData")



# conditional runtime -----------------------------------------------------

arr.job = TRUE # is this an array job or a single job
if (arr.job) {
  # get job number from array job env. var.
  job = Sys.getenv("PBS_ARRAYID")
  job = as.numeric(job)
}

# should the models be fit with species formula?
# should be one of "species", "nospecies", "both"
# should be passed from the PBS batch job
if ("both" %in% commandArgs()) {species.model = "both"}
if ("species" %in% commandArgs()) {species.model = "species"}
if ("nospecies" %in% commandArgs()) {species.model = "nospecies"}


# use all communities?
subset.data = FALSE # specify T/F to subset data
subset.size = 3000 # specifcy random subset size

# define nRCP (number of communities/RCPs)
if (arr.job) {
  #   # make the first nRCP to test = 2 (can't solve nRCP=1)
  #   nRCP = job+1
  # temp hack code to run multiple starts
  starts = rep(c(7:15),1)
  nRCP = starts[job]
} else {
  # nominal nRCP if not array job
  nRCP = 3
}



# prep model data ---------------------------------------------------------

#subset if that's happening
if (subset.data) {
  # specify the  subset to use in the analysis
  set.seed(subset.size)
  sample.sites <- sample(nrow(covariates.species), subset.size)
  covariates.species <- covariates.species[sample.sites,]
  print(paste0("Successfully subsetted [",subset.size,"] random sites"))
} else {
  print("No subsetting performed")
}


# define where abundance data starts in covariates.species, test with names(covariates.species)[n.abund]
n.abund = 17

# choose species to model with
# create a list of species with occurance > n to use for modelling
species.n = 500 # change this to desired minimum occurance
species.count = data.frame(count=sort(colSums(covariates.species[,n.abund:ncol(covariates.species)]), decreasing=T))
species.count$species = row.names(species.count)
model.species.vector = species.count$species[species.count$count>species.n]
model.species.string = paste0("cbind(", paste(model.species.vector, collapse=","),")")

# choose which covariates to use
model.covariates.vector = c("lf_rough1000_f","lf_dems1s_f","cw_rain_sumwin_f","cw_precipann_f", 
                            "ct_tempmtcp_f","defAnnual","aetAnnual","finalsurplusAnn",
                            "cw_prescott_f","fltm2_ann","SilicIndex")
model.covariates.string <- character(length(model.covariates.vector)*2)
for (i in 1:length(model.covariates.vector)){
  pos2 <- i*2
  pos1 <- pos2-1
  model.covariates.string[pos1] <- paste0(model.covariates.vector[i],".1")
  model.covariates.string[pos2] <- paste0(model.covariates.vector[i],".2")
}
model.covariates.string <- paste0(model.covariates.string, collapse = "+")

# calculate quadratic polynomial cols
covar.data = covariates.species[,model.covariates.vector]
# calculate quadratic polynomial cols
covar.data = data.frame(poly(covar.data$lf_rough1000_f, 2),
                        poly(covar.data$lf_dems1s_f, 2),
                        poly(covar.data$cw_rain_sumwin_f, 2),
                        poly(covar.data$cw_precipann_f, 2),
                        poly(covar.data$ct_tempmtcp_f, 2),
                        poly(covar.data$defAnnual, 2),
                        poly(covar.data$aetAnnual, 2),
                        poly(covar.data$finalsurplusAnn, 2),
                        poly(covar.data$cw_prescott_f, 2),
                        poly(covar.data$fltm2_ann, 2),
                        poly(covar.data$SilicIndex, 2))
names(covar.data) = c("lf_rough1000_f.1","lf_rough1000_f.2",
                      "lf_dems1s_f.1","lf_dems1s_f.2",
                      "cw_rain_sumwin_f.1","cw_rain_sumwin_f.2",
                      "cw_precipann_f.1","cw_precipann_f.2",
                      "ct_tempmtcp_f.1","ct_tempmtcp_f.2",
                      "defAnnual.1","defAnnual.2",
                      "aetAnnual.1","aetAnnual.2",
                      "finalsurplusAnn.1","finalsurplusAnn.2",
                      "cw_prescott_f.1","cw_prescott_f.2",
                      "fltm2_ann.1","fltm2_ann.2",
                      "SilicIndex.1","SilicIndex.2")


## convert categorical variables to factors?

# generate model data
model.data = data.frame(covariates.species[,model.species.vector], covar.data)

# define model form
RCP.form = paste0(model.species.string,"~","1","+",model.covariates.string)

# # add column containing factor for form.spp
# # e.g. model.data$observer = as.factor(covariates.species$Observers)
# model.data$score.method = as.factor(covariates.species$Species.score.method)
# model.data$date.int = scale(as.integer(covariates.species$Date), center=T, scale=T)
# 
# # define species form
# # e.g. "~Observer"
# species.form = "~score.method+date.int"

# record site order
site.names = covariates.species$SiteNo

# clear unused variables
rm(covariates.species, covar.data, species.count)
gc()

# conditional processing to subset data to subset - ensure sitename column is correctly specified



# fit mixture models ------------------------------------------------------
my.cont = list(maxit=3000, penalty=0.0001, penalty.tau=10, penalty.gamma=10)

if (species.model %in% c("both", "nospecies")) {
  species.form <- NULL
  tic = proc.time()
  fit.regi = regimix(form.RCP=RCP.form, form.spp=species.form, data=model.data, nRCP=nRCP, 
                     dist="Bernoulli", control=my.cont, inits="noPreClust", titbits=TRUE)
  toc = proc.time()
  
  # write model fit stats
  modelStats=list(sites=site.names, covariates=model.covariates.vector, species=model.species.vector, 
                  SppMin=species.n, SppN=fit.regi$S, nRCP=fit.regi$nRCP, runtime=round((toc-tic)[3]/60),
                  AIC=fit.regi$AIC, BIC=fit.regi$BIC, postProbs=fit.regi$postProbs, logl=fit.regi$logl, 
                  coefs=fit.regi$coefs, species.form=species.form, penalties=unlist(my.cont), conv=fit.regi$conv)
  save(modelStats, file=paste0("results/nospec/RegimixStats.n",fit.regi$n,
                               ".rcp",fit.regi$nRCP,".s",fit.regi$S,round(fit.regi$logl),".RData"))
}