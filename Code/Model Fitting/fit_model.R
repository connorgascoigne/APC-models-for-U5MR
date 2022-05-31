# packages ----

library(haven) # data import
library(latex2exp) # latex formula in plots
library(tidyverse) # data wrangling
library(SUMMER) # U5MR package
library(foreign) # importing data
library(INLA) # fitting models
library(patchwork) # putting ggplots together
library(farver) # for plots with latex labels
library(rgdal) # spatial modelling
library(ggforce) # multi page facet

# DHS  ----

## path to all Kenya folders ----

kenyaPath <- 'C:/Users/cg863/OneDrive - University of Bath/Bath PhD/Year 4/DHS/kenya2014DHS/'
# kenyaPath <- '/beegfs/scratch/user/r/cg863/code/Kenya2014DHS/'

## import functions ----

source(paste0(kenyaPath, 'Code/functions.R'))


## import data ----

# read in data
kenyaBirths <-
  haven::read_dta(paste0(kenyaPath, 'birthRecode/KEBR72DT/KEBR72FL.DTA')) %>%
  as.data.frame %>%
  mutate(b5 = haven::as_factor(b5),
         v024 = haven::as_factor(v024))

## import spatial polygons ----

### points
#### redefine regions based on clusters
kenyaPoints <- rgdal::readOGR(paste0(kenyaPath, '/subnationalBoundaries/pnts'), layer = 'KEGE71FL', verbose = FALSE)
#### poly for adjacency matrix
kenyaPoly <- rgdal::readOGR(paste0(kenyaPath, '/subnationalBoundaries/shps'), layer = 'sdr_subnational_boundaries2', verbose = FALSE)
Amat <- SUMMER::getAmat(kenyaPoly, names = kenyaPoly$REGNAME)

## import urban/rural sampling proportions ----

nationalWeights <- readRDS(paste0(kenyaPath, 'proportions/kenya_national_weights_2006_2014.rds'))
subnationalWeights <- readRDS(paste0(kenyaPath, 'proportions/kenya_subnational_weights_2006_2014.rds'))
admin1Proportions <- readRDS(paste0(kenyaPath,'proportions/kenya_admin1_proprtions_2006_2014.rds'))

# model fitting ----

## data ----

# when not using region X urbanicity strata
datTemp <- getAPCdata(data = kenyaBirths, year.cut = seq(2006, 2015, by = 1), strata = c('v025'))
dat <- region.change(data = datTemp, geoPoints = kenyaPoints)

## command line arguments ----

## First read in the arguments listed at the command line
args <- commandArgs(trailingOnly = TRUE)

## args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args) == 0){
  print("WARNING: No arguments supplied.")
  type.st <- 4
  mod <- c('apc', 'ap', 'ac')[1]
} else {
  type.st <- as.numeric(args[1]) # interaction type
  mod <- c('apc', 'ap', 'ac')[as.numeric(args[2])] # model choice from supplied index
}

# based on model choice what is slopeDrop term
if(mod == 'apc'){
  slopeDrop = 'c'
} else {
  slopeDrop = NULL
}

## create directories for results ----
### has to be done in path folder order

### All results directory

if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results'))
}

### Interaction type directory

if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Type ', type.st))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Type ', type.st))
}

### Model directory

if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Type ', type.st, '/', mod))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Type ', type.st, '/', mod))
}

## model fit ----

fit <-
  reparam.model(data = dat,
                type.st = type.st, mod = mod, slopeDrop = slopeDrop, # command line args
                inla.mode = 'experimental', is.strata = TRUE, intercept = FALSE, Amat = Amat, 
                control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE), # expensive fit
                # control.compute = list(config = TRUE), # cheap fit
                safe = TRUE)

## saving model fit ----

saveRDS(fit, file = paste0(kenyaPath, 'Code/Results/Type ', type.st, '/', mod, '/model_fit.rds'))

## saving model summary ----

sink(paste0(kenyaPath, 'Code/Results/Type ', type.st, '/', mod, '/summary.txt'))
print(summary(fit$fit))
sink()  # returns output to the console

## posterior summary ----

res <- getSamples(fit, strata.weights = NULL, save.u5m.draws = TRUE)

## saving posterior summary ----

saveRDS(res, file = paste0(kenyaPath, 'Code/Results/Type ', type.st, '/', mod, '/model_res.rds'))