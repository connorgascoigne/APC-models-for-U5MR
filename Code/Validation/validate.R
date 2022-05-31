# packages ----

library(haven) # data import
library(tidyverse) # just bare stuff for R
library(SUMMER) # U5MR package from wakeyyy
library(foreign) # importing data
library(INLA) # fitting models
library(rgdal) # spatial modelling
library(patchwork) # combining ggplots

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
  i <- 1
} else {
  type.st <- as.numeric(args[1]) # interaction type
  mod <- c('apc', 'ap', 'ac')[as.numeric(args[2])] # model choice from supplied index
  i <- as.numeric(args[3]) # region to leave out
}

# based on model choice what is slopeDrop term
if(mod == 'apc'){
  slopeDrop = 'c'
} else {
    slopeDrop = NULL
  }

## data modification ---- 

region_list <- dat$region %>% unique %>% as.vector %>% sort
last_year <- dat$period %>% unique %>% as.vector %>% max

# print in console for location check
print(paste0('interaction type: ', type.st))
print(paste0('model: ', mod))
print(paste0('region dropped: ', region_list[i]))

# change the data for each region
dataMod <-
  dat  %>%
  # complete so all region/strata combiantions have each APC combination
  tidyr::complete(nesting(region, urban, strata),
                  nesting(age, ageMid, period, perMid, cohort)) %>%
  # if the region and year combo is present replace with NA
  mutate(y = ifelse(region == region_list[i] & period == last_year, NA, y),
         N = ifelse(region == region_list[i] & period == last_year, NA, N)) %>%
  as.data.frame()

## model fit ----

fit <-
  reparam.model(data = dataMod,
                type.st = type.st, mod = mod, slopeDrop = slopeDrop, # command line args
                inla.mode = 'experimental', is.strata = TRUE, intercept = FALSE, Amat = Amat,
                safe = TRUE, control.compute = list(config = TRUE))

## posterior summaries ----

### full ----

# results
res <- getSamples(fit, strata.weights = subnationalWeights, save.u5m.draws = TRUE)

### reduced ----

# reduced version of the summary
CI <- 0.95
# credible intervals range
lowerCI <- (1 - CI)/2
upperCI <- 1 - lowerCI

# subnational weights
## includes proportion of each region for country as whole
subnationalWeights3 <-
  # join subational weights and admin1 region proportions
  subnationalWeights %>%
  # rename urban to strata
  mutate(strata = urban) %>%
  select(region, strata, period, prop)

# final 
reducedRes <- 
  # join the complete weights with u5m stratified results
  left_join(subnationalWeights3, res$u5m.draws, by = c('region', 'strata', 'period')) %>%
  # remove strata that do not exist
  filter(!(region %in% c('mombasa', 'nairobi') & strata == 'rural')) %>% 
  # label appropriately, deselect unneeded columns, organise 
  mutate(prop = prop.x) %>%
  select(-prop.x, -prop.y) %>%
  # filter out for region and year of interest 
  filter(region == region_list[i], period == last_year) %>%
  relocate(c('region', 'strata', 'period', 'prop'), .before = 'theta:1') %>%
  # multiply stratified u5m by proportion
  mutate(across(starts_with('theta:'), ~ .x * prop)) %>%
  # sum over proportions
  ungroup %>% group_by(region, period) %>%
  summarise(across(starts_with('theta:'), sum)) %>%
  # find the statistical summaries and remove individual thetas
  ungroup %>%
  mutate(select(., starts_with('theta:')) %>%
           apply(., 1, my.summary, lowerCI = lowerCI, upperCI = upperCI) %>%
           lapply(., data.frame) %>%
           do.call(rbind, .)) %>%
  mutate(width = upper - lower) %>%
  as.data.frame()

## saving posterior summaries ----

### creating directory ----

# results folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results'))
}

# validation folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation'))
}

# validation folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation'))
}

# results folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation/Results'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation/Results'))
}

# Type folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation/Results/Type ', type.st))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation/Results/Type ', type.st))
}

# Model folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation/Results/Type ', type.st, '/', mod))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation/Results/Type ', type.st, '/', mod))
}

# Full folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation/Results/Type ', type.st, '/', mod, '/Full'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation/Results/Type ', type.st, '/', mod, '/Full'))
}

# Reduced folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation/Results/Type ', type.st, '/', mod, '/Reduced'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation/Results/Type ', type.st, '/', mod, '/Reduced'))
}

saveRDS(res, file = paste0(kenyaPath, 'Code/Validation/Results/Type ', type.st, '/', mod, '/Full/leave_out_', region_list[i], '_', last_year, '.rds'))
saveRDS(reducedRes, file = paste0(kenyaPath, 'Code/Validation/Results/Type ', type.st, '/', mod, '/Reduced/leave_out_', region_list[i], '_', last_year, '.rds'))


