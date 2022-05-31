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
library(ggridges) # ridge plot

# DHS  ----

## path to all kenya folders ----

kenyaPath <- 'C:/Users/cg863/OneDrive - University of Bath/Bath PhD/Year 4/DHS/kenya2014DHS/'

## import functions ----

source(paste0(kenyaPath, 'Code/functions.R'))

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

# estimates ----

## direct ----

dir <- readRDS(file = paste0(kenyaPath, 'Code/Results/SUMMER/direct_estimate_res.rds'))

# validation results ----

type.st <- 4
region_list <- colnames(Amat) %>% sort
last_year <- 2014
type.st = 4
nsamples = 1000

## results data frame ----

resColnames <- c('region', 'period', paste0('theta:', 1:nsamples), 'mean', 'variance', 'lower', 'median', 'upper', 'width')
resTemp <- array(NA, dim = c(length(region_list), length(resColnames))) %>% as.data.frame()
colnames(resTemp) <- resColnames

apcResTemp <- resTemp
apResTemp <- resTemp
acResTemp <- resTemp

for(i in 1:length(region_list)){
  
  # i <- 1
  
  apcResTemp[i,] <- 
    readRDS(file = paste0(kenyaPath, 'Code/Results/Validation/Results/Type ', type.st, '/apc/Reduced/leave_out_', region_list[i], '_', last_year, '.rds')) %>%
    mutate(region = region %>% as.character,
           period = period %>% as.character %>% as.numeric)
  
  apResTemp[i,] <- 
    readRDS(file = paste0(kenyaPath, 'Code/Results/Validation/Results/Type ', type.st, '/ap/Reduced/leave_out_', region_list[i], '_', last_year, '.rds')) %>%
    mutate(region = region %>% as.character,
           period = period %>% as.character %>% as.numeric)
  
  acResTemp[i,] <- 
    readRDS(file = paste0(kenyaPath, 'Code/Results/Validation/Results/Type ', type.st, '/ac/Reduced/leave_out_', region_list[i], '_', last_year, '.rds')) %>%
    mutate(region = region %>% as.character,
           period = period %>% as.character %>% as.numeric)
  
}

## organising predictions ----

apcRes <-
  apcResTemp %>%
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period))

apRes <-
  apResTemp %>%
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period))

acRes <-
  acResTemp %>%
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period))

## organising 'truth' ----

dirTrue <-
  dir %>%
  filter(years == last_year,
         region != 'All') %>%
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(years),
         true = mean,
         logit.true = logit.est) %>%
  select(region, period, true, logit.true)


# scores ----

## bias and MSE from point estimate (PE) ----

### logit scale ----

logit.apcPE <-
  apcResTemp %>%
  # recode region and period
  # put U5MR-PD on logit scale 
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period),
         across(starts_with('theta:'), ~SUMMER::logit(.x))) %>%
  # drop summaries
  select(-mean, -variance, -upper, -median, -lower) %>%
  # recalc summaries on logit scale
  mutate(select(., starts_with('theta:')) %>%
           apply(., 1, my.summary, lowerCI = lowerCI, upperCI = upperCI) %>%
           lapply(., data.frame) %>%
           do.call(rbind, .)) %>%
  # join direct estimates
  left_join(., dirTrue, by = c('region', 'period')) %>% 
  select(region, period, median, variance, logit.true) %>%
  mutate(logit.bias = median - logit.true,
         logit.variance = variance,
         model = 'APC') %>%
  select(region, period, logit.bias, logit.variance, model) %>%
  # drop na regions
  drop_na() %>%
  group_by(model) %>%
  # average over all regions, square the bias and calc MSE
  summarise(logit.bias = mean(logit.bias),
            logit.variance = mean(logit.variance)) %>%
  mutate(logit.bias.squared = logit.bias^2) %>%
  mutate(logit.MSE = logit.bias.squared + logit.variance)

logit.apPE <-
  apResTemp %>%
  # recode region and period
  # put U5MR-PD on logit scale 
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period),
         across(starts_with('theta:'), ~SUMMER::logit(.x))) %>%
  # drop summaries
  select(-mean, -variance, -upper, -median, -lower) %>%
  # recalc summaries on logit scale
  mutate(select(., starts_with('theta:')) %>%
           apply(., 1, my.summary, lowerCI = lowerCI, upperCI = upperCI) %>%
           lapply(., data.frame) %>%
           do.call(rbind, .)) %>%
  # join direct estimates
  left_join(., dirTrue, by = c('region', 'period')) %>% 
  select(region, period, median, variance, logit.true) %>%
  mutate(logit.bias = median - logit.true,
         logit.variance = variance,
         model = 'AP') %>%
  select(region, period, logit.bias, logit.variance, model) %>%
  # drop na regions
  drop_na() %>%
  group_by(model) %>%
  # average over all regions, square the bias and calc MSE
  summarise(logit.bias = mean(logit.bias),
            logit.variance = mean(logit.variance)) %>%
  mutate(logit.bias.squared = logit.bias^2) %>%
  mutate(logit.MSE = logit.bias.squared + logit.variance)

logit.acPE <-
  acResTemp %>%
  # recode region and period
  # put U5MR-PD on logit scale 
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period),
         across(starts_with('theta:'), ~SUMMER::logit(.x))) %>%
  # drop summaries
  select(-mean, -variance, -upper, -median, -lower) %>%
  # recalc summaries on logit scale
  mutate(select(., starts_with('theta:')) %>%
           apply(., 1, my.summary, lowerCI = lowerCI, upperCI = upperCI) %>%
           lapply(., data.frame) %>%
           do.call(rbind, .)) %>%
  # join direct estimates
  left_join(., dirTrue, by = c('region', 'period')) %>% 
  select(region, period, median, variance, logit.true) %>%
  mutate(logit.bias = median - logit.true,
         logit.variance = variance,
         model = 'AC') %>%
  select(region, period, logit.bias, logit.variance, model) %>%
  # drop na regions
  drop_na() %>%
  group_by(model) %>%
  # average over all regions, square the bias and calc MSE
  summarise(logit.bias = mean(logit.bias),
            logit.variance = mean(logit.variance)) %>%
  mutate(logit.bias.squared = logit.bias^2) %>%
  mutate(logit.MSE = logit.bias.squared + logit.variance)

### expit scale ----

expit.apcPE <-
  apcResTemp %>%
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period)) %>%
  left_join(., dirTrue, by = c('region', 'period')) %>% 
  # on expit scale, bias against the DE and variance recalced
  mutate(expit.bias = median - true,
         expit.variance = variance,
         model = 'APC') %>%
  select(region, period, mean, median, variance, true, expit.bias, expit.variance, model) %>%
  # drop na regions
  drop_na() %>%
  group_by(model) %>%
  # average over all regions, square the bias and calc MSE
  summarise(expit.bias = mean(expit.bias),
            expit.variance = mean(expit.variance)) %>%
  mutate(expit.bias.squared = expit.bias^2) %>%
  mutate(expit.MSE = expit.bias.squared + expit.variance)

expit.apPE <-
  apResTemp %>%
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period)) %>%
  left_join(., dirTrue, by = c('region', 'period')) %>%  
  # on expit scale, bias against the DE and variance recalced
  mutate(expit.bias = median - true,
         expit.variance = variance,
         model = 'AP') %>%
  select(region, period, mean, median, variance, true, expit.bias, expit.variance, model) %>%
  drop_na() %>%
  # drop na regions
  group_by(model) %>%
  # average over all regions, square the bias and calc MSE
  summarise(expit.bias = mean(expit.bias),
            expit.variance = mean(expit.variance)) %>%
  mutate(expit.bias.squared = expit.bias^2) %>%
  mutate(expit.MSE = expit.bias.squared + expit.variance)

expit.acPE <-
  acResTemp %>%
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period)) %>%
  left_join(., dirTrue, by = c('region', 'period')) %>% 
  # on expit scale, bias against the DE and variance recalced
  mutate(expit.bias = median - true,
         expit.variance = variance,
         model = 'AC') %>%
  select(region, period, mean, median, variance, true, expit.bias, expit.variance, model) %>%
  # drop na regions
  drop_na() %>%
  group_by(model) %>%
  # average over all regions, square the bias and calc MSE
  summarise(expit.bias = mean(expit.bias),
            expit.variance = mean(expit.variance)) %>%
  mutate(expit.bias.squared = expit.bias^2) %>% 
  mutate(expit.MSE = expit.bias.squared + expit.variance)

### logit and expit together ----

expit.peScores <- 
  rbind(expit.apcPE, expit.apPE, expit.acPE) %>%
  select(expit.bias.squared, expit.variance, expit.MSE) %>%
  mutate(expit.bias.squared = expit.bias.squared * 1e5,
         expit.variance = expit.variance * 1e5,
         expit.MSE = expit.MSE * 1e5)

logit.peScores <- 
  rbind(logit.apcPE, logit.apPE, logit.acPE) %>%
  select(logit.bias.squared, logit.variance, logit.MSE) %>%
  mutate(logit.bias.squared = logit.bias.squared * 1e3,
         logit.variance = logit.variance * 1e3,
         logit.MSE = logit.MSE * 1e3)

names <- data.frame(model = c('APC', 'AP', 'AC'))

peScores <- cbind(names, expit.peScores, logit.peScores)

print(xtable::xtable(peScores, type = 'latex'), include.rownames = FALSE, digits = 2)

## cpo ----

# direct estimate

dirTemp <-
  dir %>%
  mutate(period = as.numeric(years),
         directLogitEstimate = logit.est,
         directVariance = var.est) %>%
  select(region, period, directVariance)

# apc 

apcCPOTemp <-
  apcResTemp %>%
  mutate(predMean = mean) %>%
  select(region, period, predMean, starts_with('theta:')) %>%
  # join direct estimate
  left_join(., dirTemp, by = c('region', 'period')) %>%
  # theta here are place holders
  select(region, period, predMean, directVariance, starts_with('theta:')) %>%
  # replace the NAs with 0
  mutate(directVariance = replace_na(directVariance, 0),
         across(starts_with('theta:'),~ NA))

set.seed(123)

for(i in 1:length(region_list)){
  
  # create a sample from the asymptotic distribution
  apcCPOTemp[i,5:ncol(apcCPOTemp)] <- 
    rnorm(n = nsamples, 
          mean = SUMMER::logit(apcCPOTemp$predMean[i]), 
          sd = sqrt(apcCPOTemp$directVariance[i]))
  
}

apcCPO <-
  apcCPOTemp %>%
  # take the inverse link function
  mutate(across(starts_with('theta:'), ~ SUMMER::expit(.x))) %>%
  # find the point estimates (mean and quantiles)
  mutate(select(., starts_with('theta:')) %>%
           apply(., 1, my.summary, lowerCI = lowerCI, upperCI = upperCI) %>%
           lapply(., data.frame) %>%
           do.call(rbind, .)) %>%
  select(region, period, mean) %>%
  mutate(logRegionalCPO = log(mean),
         model = 'APC')

# ap

apCPOTemp <-
  apResTemp %>%
  mutate(predMean = mean) %>%
  select(region, period, predMean, starts_with('theta:')) %>%
  # join direct estimate
  left_join(., dirTemp, by = c('region', 'period')) %>%
  # theta here are place holders
  select(region, period, predMean, directVariance, starts_with('theta:')) %>%
  # replace the NAs with 0
  mutate(directVariance = replace_na(directVariance, 0),
         across(starts_with('theta:'),~ NA))

set.seed(123)

for(i in 1:length(region_list)){
  
  apCPOTemp[i,5:ncol(apCPOTemp)] <- 
    rnorm(n = nsamples, 
          mean =  SUMMER::logit(apCPOTemp$predMean[i]), 
          sd = sqrt(apCPOTemp$directVariance[i]))
  
}

apCPO <-
  apCPOTemp %>%
  mutate(across(starts_with('theta:'), ~ SUMMER::expit(.x))) %>%
  mutate(select(., starts_with('theta:')) %>%
           apply(., 1, my.summary, lowerCI = lowerCI, upperCI = upperCI) %>%
           lapply(., data.frame) %>%
           do.call(rbind, .)) %>%
  select(region, period, mean) %>%
  mutate(logRegionalCPO = log(mean),
         model = 'AP')

# ac

set.seed(123)

acCPOTemp <-
  acResTemp %>%
  mutate(predMean = mean) %>%
  select(region, period, predMean, starts_with('theta:')) %>%
  # join direct estimate
  left_join(., dirTemp, by = c('region', 'period')) %>%
  # theta here are place holders
  select(region, period, predMean, directVariance, starts_with('theta:')) %>%
  # replace the NAs with 0
  mutate(directVariance = replace_na(directVariance, 0),
         across(starts_with('theta:'),~ NA))


for(i in 1:length(region_list)){
  
  acCPOTemp[i,5:ncol(acCPOTemp)] <- 
    rnorm(n = nsamples, 
          mean =  SUMMER::logit(acCPOTemp$predMean[i]), 
          sd = sqrt(acCPOTemp$directVariance[i]))
  
}

acCPO <-
  acCPOTemp %>%
  mutate(across(starts_with('theta:'), ~ SUMMER::expit(.x))) %>%
  mutate(select(., starts_with('theta:')) %>%
           apply(., 1, my.summary, lowerCI = lowerCI, upperCI = upperCI) %>%
           lapply(., data.frame) %>%
           do.call(rbind, .)) %>%
  select(region, period, mean) %>%
  mutate(logRegionalCPO = log(mean),
         model = 'AC')

# scores 

# find the overall score by summing the log of the CPOS
allCPO <- 
  rbind(apcCPO, apCPO, acCPO) %>%
  group_by(model) %>%
  summarise(LCPO = sum(logRegionalCPO)) %>%
  as.data.frame()




## coverage ----
# does my estimate lie within the 95% CI for direct estimate

dirTemp <-
  dir %>%
  mutate(period = as.numeric(years)) %>%
  select(region, period, mean, lower, upper)


apcCoverage <-
  apcResTemp %>%
  select(region, period, starts_with('theta:')) %>%
  left_join(., dirTemp, by = c('region', 'period')) %>%
  mutate(across(starts_with('theta:'), ~ if_else( lower < .x & .x < upper, 1, 0))) %>%
  mutate(coverage = rowMeans(select(., starts_with('theta:')), na.rm = TRUE),
         model = 'APC') %>%
  relocate(c('model', 'mean', 'lower', 'upper', 'coverage'), .before = 'theta:1') %>%
  na.omit() %>%
  group_by(model) %>%
  summarise(coverage = mean(coverage))

apCoverage <-
  apResTemp %>%
  select(region, period, starts_with('theta:')) %>%
  left_join(., dirTemp, by = c('region', 'period')) %>%
  mutate(across(starts_with('theta:'), ~ if_else( lower < .x & .x < upper, 1, 0))) %>%
  mutate(coverage = rowMeans(select(., starts_with('theta:')), na.rm = TRUE),
         model = 'AP') %>%
  relocate(c('model', 'mean', 'lower', 'upper', 'coverage'), .before = 'theta:1') %>%
  na.omit() %>%
  group_by(model) %>%
  summarise(coverage = mean(coverage))

acCoverage <-
  acResTemp %>%
  select(region, period, starts_with('theta:')) %>%
  left_join(., dirTemp, by = c('region', 'period')) %>%
  mutate(across(starts_with('theta:'), ~ if_else( lower < .x & .x < upper, 1, 0))) %>%
  mutate(coverage = rowMeans(select(., starts_with('theta:')), na.rm = TRUE),
         model = 'AC') %>%
  relocate(c('model', 'mean', 'lower', 'upper', 'coverage'), .before = 'theta:1') %>%
  na.omit() %>%
  group_by(model) %>%
  summarise(coverage = mean(coverage))

allCoverage <- 
  rbind(apcCoverage, apCoverage, acCoverage) %>%
  mutate(coverage = coverage * 100)

print(xtable::xtable(allCoverage, type = 'latex'), include.rownames = FALSE, digits = 2)

## combining all score ----

finalResults <- 
  left_join(peScores, allCoverage, by = 'model') %>%
  left_join(., allCPO, by = 'model')

print(xtable::xtable(t(finalResults %>% select(-model)), type = 'latex'), digits = 2)

# plots ----

## directory ----

# results folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results'))
}

# validation folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation'))
}

# plots folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation/Plots'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation/Plots'))
}

# type folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st))
}

# predicition ridgeplot folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Validation Ridge Plots'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Validation Ridge Plots'))
}

# cpo ridgeplot folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/CPO Ridge Plots'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/CPO Ridge Plots'))
}

# Subnational Maps
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Subnational Maps'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Subnational Maps'))
}

## U5MR-PD ridge histogram ----

apcRidgeRes <- 
  apcRes %>%
  select(region, starts_with('theta:')) %>%
  pivot_longer(cols = starts_with('theta:'),
               names_to = 'theta',
               values_to = 'estimate')

apRidgeRes <- 
  apRes %>%
  select(region, starts_with('theta:')) %>%
  pivot_longer(cols = starts_with('theta:'),
               names_to = 'theta',
               values_to = 'estimate')

acRidgeRes <- 
  acRes %>%
  select(region, starts_with('theta:')) %>%
  pivot_longer(cols = starts_with('theta:'),
               names_to = 'theta',
               values_to = 'estimate')

apcRidge<-
  ggplot(apcRidgeRes, aes(x = estimate * 1000, y = str_to_title(region), fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  ggplot2::scale_fill_viridis_c(name = 'Predicted U5MR \n(per 1000 live births)', direction = -1, option = 'D') +
  ggplot2::scale_y_discrete(limits = rev) +
  ggplot2::labs(y = 'Region') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  my.theme(legend.position = 'bottom',
           axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank())

apRidge<-
  ggplot(apRidgeRes, aes(x = estimate * 1000, y = str_to_title(region), fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  ggplot2::scale_fill_viridis_c(name = 'Predicted U5MR \n(per 1000 live births)', direction = -1, option = 'D') +
  ggplot2::scale_y_discrete(limits = rev) +
  ggplot2::labs(y = 'Region') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  my.theme(legend.position = 'bottom',
           axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank())

acRidge<-
  ggplot(acRidgeRes, aes(x = estimate * 1000, y = str_to_title(region), fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  ggplot2::scale_fill_viridis_c(name = 'Predicted U5MR \n(per 1000 live births)', direction = -1, option = 'D') +
  ggplot2::scale_y_discrete(limits = rev) +
  ggplot2::labs(y = 'Region') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  my.theme(legend.position = 'bottom',
           axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank())

ggsave(apcRidge, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Validation Ridge Plots/apcValidationRidgePlotType', type.st,  '.png'), height = 10, width = 10)
ggsave(apRidge, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Validation Ridge Plots/apValidationRidgePlotType', type.st,  '.png'), height = 10, width = 10)
ggsave(acRidge, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Validation Ridge Plots/acValidationRidgePlotType', type.st,  '.png'), height = 10, width = 10)

## CPO ridge histogram ----

apcCPOridgeRes <- 
  apcCPOTemp %>%
  # take the inverse link function
  mutate(across(starts_with('theta:'), ~ SUMMER::expit(.x))) %>%
  select(region, starts_with('theta:')) %>%
  pivot_longer(cols = starts_with('theta:'),
               names_to = 'theta',
               values_to = 'estimate')

apCPOridgeRes <- 
  apCPOTemp %>%
  # take the inverse link function
  mutate(across(starts_with('theta:'), ~ SUMMER::expit(.x))) %>%
  select(region, starts_with('theta:')) %>%
  pivot_longer(cols = starts_with('theta:'),
               names_to = 'theta',
               values_to = 'estimate')

acCPOridgeRes <- 
  acCPOTemp %>%
  # take the inverse link function
  mutate(across(starts_with('theta:'), ~ SUMMER::expit(.x))) %>%
  select(region, starts_with('theta:')) %>%
  pivot_longer(cols = starts_with('theta:'),
               names_to = 'theta',
               values_to = 'estimate')

apcCPOridge<-
  ggplot(apcCPOridgeRes, aes(x = estimate * 1000, y = str_to_title(region), fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  ggplot2::scale_fill_viridis_c(name = 'Predicted U5MR \n(per 1000 live births)', direction = -1, option = 'D') +
  ggplot2::scale_y_discrete(limits = rev) +
  ggplot2::labs(y = 'Region') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  ggplot2::geom_segment(data = dirTrue, aes(x = 1000*true, xend = 1000*true, y = as.numeric(region), yend = as.numeric(region) + .9), color = 'red') +
  my.theme(legend.position = 'bottom',
           axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank())

apCPOridge<-
  ggplot(apCPOridgeRes, aes(x = estimate * 1000, y = str_to_title(region), fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  ggplot2::scale_fill_viridis_c(name = 'Predicted U5MR \n(per 1000 live births)', direction = -1, option = 'D') +
  ggplot2::scale_y_discrete(limits = rev) +
  ggplot2::labs(y = 'Region') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  ggplot2::geom_segment(data = dirTrue, aes(x = 1000*true, xend = 1000*true, y = as.numeric(region), yend = as.numeric(region) + .9), color = 'red') +
  my.theme(legend.position = 'bottom',
           axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank())

acCPOridge<-
  ggplot(acCPOridgeRes, aes(x = estimate * 1000, y = str_to_title(region), fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  ggplot2::scale_fill_viridis_c(name = 'Predicted U5MR \n(per 1000 live births)', direction = -1, option = 'D') +
  ggplot2::scale_y_discrete(limits = rev) +
  ggplot2::labs(y = 'Region') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  ggplot2::geom_segment(data = dirTrue, aes(x = 1000*true, xend = 1000*true, y = as.numeric(region), yend = as.numeric(region) + .9), color = 'red') +
  my.theme(legend.position = 'bottom',
           axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank())

# ggsave(apcCPORidge, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/CPO Ridge Plots/apcCPORidgePlotType', type.st,  '.png'), height = 10, width = 10)
# ggsave(apCPORidge, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/CPO Ridge Plots/apCPORidgePlotType', type.st,  '.png'), height = 10, width = 10)
# ggsave(acCPORidge, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/CPO Ridge Plots/acCPORidgePlotType', type.st,  '.png'), height = 10, width = 10)

## map plots ----

apcMapRes <- 
  apcRes %>%
  select(-starts_with('theta:')) %>%
  mutate(region = region %>% as.character) %>%
  as.data.frame

apMapRes <- 
  apRes %>%
  select(-starts_with('theta:')) %>%
  mutate(region = region %>% as.character) %>%
  as.data.frame

acMapRes <- 
  acRes %>%
  select(-starts_with('theta:')) %>%
  mutate(region = region %>% as.character) %>%
  as.data.frame

# apc
apcMedianMap <-
  my.map.plot(data = apcMapRes,
              geo = kenyaPoly,
              variables = 'period',
              value = 'median',
              by.data = 'region',
              by.geo = 'REGNAME',
              ncol = 3,
              per1000 = TRUE,
              direction = - 1,
              color.option = 'D',
              legend.label = 'Predicted U5MR \n(per 1000 live births)') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  my.map.theme(legend.position = 'bottom', 
               strip.text.x = element_blank(),
               legend.title.align = 0.5)


apcWidthMap <-
  my.map.plot(data = apcMapRes,
              geo = kenyaPoly,
              variables = 'period',
              value = 'width',
              by.data = 'region',
              by.geo = 'REGNAME',
              ncol = 3,
              per1000 = TRUE,
              direction = - 1,
              color.option = 'A',
              legend.label = 'Width of 95% CI') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  my.map.theme(legend.position = 'bottom', 
               strip.text.x = element_blank(),
               legend.title.align = 0.5)

# ap
apMedianMap <-
  my.map.plot(data = apMapRes,
              geo = kenyaPoly,
              variables = 'period',
              value = 'median',
              by.data = 'region',
              by.geo = 'REGNAME',
              ncol = 3,
              per1000 = TRUE,
              direction = - 1,
              color.option = 'D',
              legend.label = 'Predicted U5MR \n(per 1000 live births)') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  my.map.theme(legend.position = 'bottom', 
               strip.text.x = element_blank(),
               legend.title.align = 0.5)


apWidthMap <-
  my.map.plot(data = apMapRes,
              geo = kenyaPoly,
              variables = 'period',
              value = 'width',
              by.data = 'region',
              by.geo = 'REGNAME',
              ncol = 3,
              per1000 = TRUE,
              direction = - 1,
              color.option = 'A',
              legend.label = 'Width of 95% CI') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  my.map.theme(legend.position = 'bottom', 
               strip.text.x = element_blank(),
               legend.title.align = 0.5)

# ac
acMedianMap <-
  my.map.plot(data = acMapRes,
              geo = kenyaPoly,
              variables = 'period',
              value = 'median',
              by.data = 'region',
              by.geo = 'REGNAME',
              ncol = 3,
              per1000 = TRUE,
              direction = - 1,
              color.option = 'D',
              legend.label = 'Predicted U5MR \n(per 1000 live births)') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  my.map.theme(legend.position = 'bottom', 
               strip.text.x = element_blank(),
               legend.title.align = 0.5)


acWidthMap <-
  my.map.plot(data = acMapRes,
              geo = kenyaPoly,
              variables = 'period',
              value = 'width',
              by.data = 'region',
              by.geo = 'REGNAME',
              ncol = 3,
              per1000 = TRUE,
              direction = - 1,
              color.option = 'A',
              legend.label = 'Width of 95% CI') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  my.map.theme(legend.position = 'bottom', 
               strip.text.x = element_blank(),
               legend.title.align = 0.5)

# apcMedianMap + apcWidthMap
# apMedianMap + apWidthMap
# acMedianMap + acWidthMap

ggsave(apcMedianMap, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Subnational Maps/apcValidationMedianMapType', type.st,  '.png'), height = 7, width = 7)
ggsave(apcWidthMap, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Subnational Maps/apcValidationWidthMapType', type.st,  '.png'), height = 7, width = 7)

ggsave(apMedianMap, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Subnational Maps/apValidationMedianMapType', type.st,  '.png'), height = 7, width = 7)
ggsave(apWidthMap, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Subnational Maps/apValidationWidthMapType', type.st,  '.png'), height = 7, width = 7)

ggsave(acMedianMap, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Subnational Maps/acValidationMedianMapType', type.st,  '.png'), height = 7, width = 7)
ggsave(acWidthMap, file = paste0(kenyaPath, 'Code/Results/Validation/Plots/Type ', type.st, '/Subnational Maps/acValidationWidthMapType', type.st,  '.png'), height = 7, width = 7)


