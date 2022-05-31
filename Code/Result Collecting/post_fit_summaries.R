# packages ----

library(haven) # data import
library(tidyverse) # just bare stuff for R
library(SUMMER) # U5MR package from wakeyyy
library(foreign) # importing data
library(INLA) # fitting models
library(rgdal) # spatial modelling
library(patchwork) # combining ggplots

# DHS  ----

## path to all kenya folders ----

kenyaPath <- 'C:/Users/cg863/OneDrive - University of Bath/Bath PhD/Year 4/DHS/kenya2014DHS/'

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

## APC ----

### APC data ----

# when not using region X urbanicity strata
datTemp <- getAPCdata(data = kenyaBirths, year.cut = seq(2006, 2015, by = 1), strata = c('v025'))
dat <- region.change(data = datTemp, geoPoints = kenyaPoints)

### APC models ----

options(scipen = 100000)
# options(scipen = 0) # reset
digits = 3

#### Type I ----

load(file = paste0(kenyaPath, 'Code/Model Fits/Type_I_fits.rda'))

apcSummary <- summary(apcFit$fit, digits = digits)
apSummary <- summary(apFit$fit, digits = digits)
acSummary <- summary(acFit$fit, digits = digits)

sink(paste0(kenyaPath, 'Code/Model Fits/Type_I_summaries.txt'))
print('APC')
print(apcSummary)
print('AP')
print(apSummary)
print('AC')
print(acSummary)
sink()  # returns output to the console

rm(apcFit, apcRes, apFit, apRes, acFit, acRes)

#### Type II ----

load(file = paste0(kenyaPath, 'Code/Model Fits/Type_II_fits.rda'))

apcSummary <- summary(apcFit$fit, digits = digits)
apSummary <- summary(apFit$fit, digits = digits)
acSummary <- summary(acFit$fit, digits = digits)

sink(paste0(kenyaPath, 'Code/Model Fits/Type_II_summaries.txt'))
print('APC')
print(apcSummary)
print('AP')
print(apSummary)
print('AC')
print(acSummary)
sink()  # returns output to the console

rm(apcFit, apcRes, apFit, apRes, acFit, acRes)

#### Type III ----

load(file = paste0(kenyaPath, 'Code/Model Fits/Type_III_fits.rda'))

apcSummary <- summary(apcFit$fit, digits = digits)
apSummary <- summary(apFit$fit, digits = digits)
acSummary <- summary(acFit$fit, digits = digits)

sink(paste0(kenyaPath, 'Code/Model Fits/Type_III_summaries.txt'))
print('APC')
print(apcSummary)
print('AP')
print(apSummary)
print('AC')
print(acSummary)
sink()  # returns output to the console

rm(apcFit, apcRes, apFit, apRes, acFit, acRes)

#### Type IV ----

load(file = paste0(kenyaPath, 'Code/Model Fits/Type_IV_fits.rda'))

apcSummary <- summary(apcFit$fit, digits = digits)
apSummary <- summary(apFit$fit, digits = digits)
acSummary <- summary(acFit$fit, digits = digits)

sink(paste0(kenyaPath, 'Code/Model Fits/Type_IV_summaries.txt'))
print('APC')
print(apcSummary)
print('AP')
print(apSummary)
print('AC')
print(acSummary)
sink()  # returns output to the console

rm(apcFit, apcRes, apFit, apRes, acFit, acRes)

# precisions as sd
# apc
1/c(29.013, 150.467, 567.367) # a
1/c(14.912, 165.977, 1968.495) # p
1/c(59.853, 784.568, 19873.478) # c
1/c(4.829, 7.923, 12.924) # r
1/c(53.716, 478.488, 8747.622) # s-t

# ap
1/c(234.970, 154.38, 69.54) #a
1/c(85.845, 1268.967, 105450.308) # p
1/c(4.876, 7.997, 13.099) # r
1/c(118.408, 616.269, 9347.391) # s-t

# ac
1/c(27.179, 146.641, 530.378) # a
1/c(95.271, 4928.381, 210282.624) # c
1/c(4.836, 7.723, 12.441) # r
1/c(426.661, 4177.737, 320582.844) # s-t