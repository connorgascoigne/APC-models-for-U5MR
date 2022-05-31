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

## path to all kenya folders ----

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

## dir ----

### data ----

smthDirDatTemp <-
  SUMMER::getBirths(data = kenyaBirths,
                    strata = 'v025',
                    surveyyear = 2014,
                    year.cut = seq(2006, 2015, by = 1)) %>%
  mutate(v005 = v005/1e6) %>%
  select(v001, v002, v005, v024, strata, time, age, died)
colnames(smthDirDatTemp) <- c('cluster', 'house', 'weight', 'region', 'strata', 'year', 'age', 'died')

smthDirDat <- region.change(data = smthDirDatTemp, geoPoints = kenyaPoints)

### estimates ----

dir <- 
  SUMMER::getDirect(births = smthDirDat,
                    years = levels(smthDirDat$year),
                    regionVar = 'region',
                    timeVar = 'year',
                    clusterVar =  '~cluster',
                    ageVar = 'age',
                    weightsVar = 'weight')

## smooth direct ----

### model fit ----

smthDirFit <- 
  smoothDirect(data = dir, 
               Amat = NULL, 
               year_label = levels(smthDirDat$year),
               year_range = c(2006, 2014), 
               type.st = 4,
               options = list(dic = TRUE, cpo = TRUE, waic = TRUE, openmp.strategy = 'default'),
               m = 1)

### posterior summary ----

smthDirRes <- getSmoothed(smthDirFit)

## smooth cluster ----

### data ----

clusterDatTemp <-
  SUMMER::getBirths(data = kenyaBirths,
                    strata = 'v025',
                    surveyyear = 2014,
                    year.cut = seq(2006, 2015, by = 1),
                    compact = T) %>%
  mutate(strata = haven::as_factor(strata)) %>%
  # more informative names
  rename(cluster = v001, region = v024, years = time, Y = died) %>%
  # select only the columns needed and arrange order
  select(cluster, region, strata, age, years, Y, total) %>%
  arrange(cluster, region, strata, age, years)

clusterDat <- region.change(clusterDatTemp, geoPoints = kenyaPoints)

### model fit ----

smthClusFit <- 
  smoothCluster(data = clusterDat,
                # stratified-national level model
                Amat = NULL, 
                family = 'betabinomial',
                year_label = levels(clusterDat$years) %>% as.numeric(),
                age.groups = levels(clusterDat$age),
                # want only one rw2 through age
                age.rw.group = rep(1, times = levels(clusterDat$age) %>% length))

### posterior summary ----

# want to aggregate over strata
## adapted form yunhans code <- this producing the same dataframe 
## confident this will be accepted into the getSmoothed function
smthClustNationalWeights <-
  nationalWeights %>%
  # change period to numeric years
  rename(years = period) %>%
  mutate(years = years %>% as.character %>% as.numeric) %>%
  # pivot to the correct formula
  tidyr::pivot_wider(names_from = urban, values_from = prop) %>%
  as.data.frame()

smthClusRes <- getSmoothed(smthClusFit, weight.strata = smthClustNationalWeights, save.draws = TRUE)

# saving ----

## create directories ----

# results folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results'))
}

# SUMMER results folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/SUMMER'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/SUMMER'))
}

## saving ----

# direct estimates
saveRDS(dir, file = paste0(kenyaPath, 'Code/Results/SUMMER/direct_estimate_res.rds'))

# smooth direct
saveRDS(smthDirFit, file = paste0(kenyaPath, 'Code/Results/SUMMER/smooth_direct_estimate_fit.rds'))
saveRDS(smthDirRes, file = paste0(kenyaPath, 'Code/Results/SUMMER/smooth_direct_estimate_res.rds'))

# smooth cluster
saveRDS(smthClusFit, file = paste0(kenyaPath, 'Code/Results/SUMMER/smooth_cluster_fit.rds'))
saveRDS(smthClusRes, file = paste0(kenyaPath, 'Code/Results/SUMMER/smooth_cluster_res.rds'))

