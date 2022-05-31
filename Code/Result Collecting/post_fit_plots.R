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

# result collecting ----

## APC ----


### models fits and results ----

## First read in the arguments listed at the command line
args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0){
  print("WARNING: No arguments supplied.")
  type.st <- 4
} else {
  type.st <- as.numeric(args[1]) # interaction type
}


apcFit <- readRDS(file = paste0(kenyaPath, 'Code/Results/Type ', type.st, '/apc/model_fit.rds'))
apFit <- readRDS(file = paste0(kenyaPath, 'Code/Results/Type ', type.st, '/ap/model_fit.rds'))
acFit <- readRDS(file = paste0(kenyaPath, 'Code/Results/Type ', type.st, '/ac/model_fit.rds'))

apcRes <- readRDS(file = paste0(kenyaPath, 'Code/Results/Type ', type.st, '/apc/model_res.rds'))
apRes <- readRDS(file = paste0(kenyaPath, 'Code/Results/Type ', type.st, '/ap/model_res.rds'))
acRes <- readRDS(file = paste0(kenyaPath, 'Code/Results/Type ', type.st, '/ac/model_res.rds'))

## SUMMER ----

### model results ----

dir <- readRDS(file = paste0(kenyaPath, 'Code/Results/SUMMER/direct_estimate_res.rds'))
smthDirRes <- readRDS(file = paste0(kenyaPath, 'Code/Results/SUMMER/smooth_direct_estimate_res.rds'))
smthClusRes <- readRDS(file = paste0(kenyaPath, 'Code/Results/SUMMER/smooth_cluster_res.rds'))

# main directories ----

# results folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results'))
}

# plots folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Plots'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Plots'))
}

# interaction type folder
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st))
}


# national results ----

## APC proportional aggregation (national) ----

CI <- 0.95
# credible intervals range
lowerCI <- (1 - CI)/2
upperCI <- 1 - lowerCI

# subnational weights
## includes proportion of each region for country as whole
subnationalWeights2 <-
  # join subational weights and admin1 region proportions
  left_join(subnationalWeights, admin1Proportions, by = c('region', 'period')) %>%
  # multiple proportions together for a urban/rural X region prop
  mutate(prop = prop.x * prop.y,
         strata = urban) %>%
  select(region, strata, period, prop)

# # check
# ## should all be 1 for each year
# subnationalWeights2 %>% group_by(period) %>% summarise(total = sum(prop))

# apc
apcNatFinal <- 
  # join the complete weights with u5m stratified results
  left_join(subnationalWeights2, apcRes$u5m.draws, by = c('region', 'strata', 'period')) %>%
  # remove strata that do not exist
  filter(!(region %in% c('mombasa', 'nairobi') & strata == 'rural')) %>% 
  # label appropriately, deselect unneeded columns, organise 
  mutate(prop = prop.x,
         region = 'All' %>% as.factor()) %>%
  select(-prop.x, -prop.y) %>%
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
  select(-starts_with('theta:'))  %>%
  mutate(type = haven::as_factor('APC')) %>%
  select(region, period, median, lower, upper, type) %>%
  as.data.frame()
# ap
apNatFinal <- 
  # join the complete weights with u5m stratified results
  left_join(subnationalWeights2, apRes$u5m.draws, by = c('region', 'strata', 'period')) %>%
  # remove strata that do not exist
  filter(!(region %in% c('mombasa', 'nairobi') & strata == 'rural')) %>% 
  # label appropriately, deselect unneeded columns, organise 
  mutate(prop = prop.x,
         region = 'All' %>% as.factor()) %>%
  select(-prop.x, -prop.y) %>%
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
  select(-starts_with('theta:'))  %>%
  mutate(type = haven::as_factor('AP')) %>%
  select(region, period, median, lower, upper, type) %>%
  as.data.frame()
# ac
acNatFinal <- 
  # join the complete weights with u5m stratified results
  left_join(subnationalWeights2, acRes$u5m.draws, by = c('region', 'strata', 'period')) %>%
  # remove strata that do not exist
  filter(!(region %in% c('mombasa', 'nairobi') & strata == 'rural')) %>% 
  # label appropriately, deselect unneeded columns, organise 
  mutate(prop = prop.x,
         region = 'All' %>% as.factor()) %>%
  select(-prop.x, -prop.y) %>%
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
  select(-starts_with('theta:'))  %>%
  mutate(type = haven::as_factor('AC')) %>%
  select(region, period, median, lower, upper, type) %>%
  as.data.frame()

## line plots ----

natRes <-
  rbind(dir %>%
          filter(region == 'All') %>%
          mutate(period = years,
                 type = haven::as_factor('Direct'),
                 median = mean) %>%
          select(region, period, median, lower, upper, type),
        smthDirRes %>%
          mutate(period = years,
                 type = haven::as_factor('Smooth Direct')) %>%
          select(region, period, median, lower, upper, type),
        # smthClusRes$overall %>%
        #   mutate(period = years,
        #          type = haven::as_factor('Age Time')) %>%
        #   select(region, period, median, lower, upper, type),
        apcNatFinal,
        apNatFinal,
        acNatFinal)

#  colour
natPlotError <-
  ggplot(natRes %>% mutate(median = 1000*median), 
         aes(x = period, y = median, color = type, group = interaction(region, type))) +
  ggplot2::geom_point(position = position_dodge(0.5)) +
  ggplot2::geom_errorbar(aes(ymin = 1000*lower, ymax = 1000*upper), width = .2, position = position_dodge(0.5)) +
  ggplot2::labs(x = 'Period', y = 'U5MR (deaths per 1000 live births)') +
  my.theme(text = element_text(size = 15),
           axis.text.x = element_text(angle = 45, hjust = 1),
           strip.text.y = element_blank(),
           legend.position = 'top',
           legend.title = element_blank())
natPlotLine <-
  ggplot(natRes %>% mutate(median = 1000*median), 
         aes(x = period, y = median, color = type, group = interaction(region, type))) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::labs(x = 'Period', y = 'U5MR (deaths per 1000 live births)') +
  my.theme(text = element_text(size = 15),
           axis.text.x = element_text(angle = 45, hjust = 1),
           strip.text.y = element_blank(),
           legend.position = 'top',
           legend.title = element_blank())

# save
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/National Line Plot'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/National Line Plot'))
}

ggsave(plot = natPlotError,
       filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/National Line Plot/natPlotErrorType', type.st, '.png'),
       height = 7, width = 7)

ggsave(plot = natPlotLine,
       filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/National Line Plot/natPlotLineType', type.st, '.png'),
       height = 7, width = 7)


## APC Estimates vs Other Estimates (national) ----

dirTemp <-
  dir %>% 
  mutate(period = years)

# data organisation
natRes <- 
  rbind(apcNatFinal, apNatFinal, acNatFinal) %>%
  left_join(., dirTemp, by = c('region', 'period')) %>%
  mutate(modEst = SUMMER::logit(median),
         dirEst = SUMMER::logit(mean)) %>%
  select(region, period, modEst, dirEst, type)

# vs direct
apcVsdirNat<-
  ggplot(natRes %>% filter(type == 'APC'), aes(y = modEst, x = dirEst)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(-3.1, -2.5) + ylim(-3.1, -2.5) +
  labs(y = 'APC Estimate (logit)', x = 'Direct Estimate (logit)') +
  my.theme(text = element_text(size = 15))

# vs ap
apcVSapNat<-
  ggplot() +
  geom_point(aes(y = SUMMER::logit(apcNatFinal$median), x = SUMMER::logit(apNatFinal$median)), size = 3) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(-3.1, -2.5) + ylim(-3.1, -2.5) +
  labs(y = 'APC Estimate (logit)', x = 'AP Estimate (logit)') +
  my.theme(text = element_text(size = 15))

# vs ac
apcVSacNat<-
  ggplot() +
  geom_point(aes(y = SUMMER::logit(apcNatFinal$median), x = SUMMER::logit(acNatFinal$median)), size = 3) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(-3.1, -2.5) + ylim(-3.1, -2.5) +
  labs(y = 'APC Estimate (logit)', x = 'AC Estimate (logit)') +
  my.theme(text = element_text(size = 15))


# save
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/National APC Vs Other Estimates'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/National APC Vs Other Estimates'))
}

# saving
height = width = 7

## apc
ggsave(plot = apcVsdirNat, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/National APC Vs Other Estimates/apcVsdirNatType', type.st, '.png'),
       height = height, width = width)

## ap
ggsave(plot = apcVSapNat, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/National APC Vs Other Estimates/apcVSapNatType', type.st, '.png'),
       height = height, width = width)
## ac
ggsave(plot = apcVSacNat, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/National APC Vs Other Estimates/apcVSacNatType', type.st, '.png'),
       height = height, width = width)




# subnational results ----

## APC proportional aggregation (subnational) ----

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

# # check
# ## should all be 1 for each year-region
# subnationalWeights %>% group_by(region, period) %>% summarise(total = sum(prop))

# apc
apcSubnatFinal <- 
  # join the complete weights with u5m stratified results
  left_join(subnationalWeights3, apcRes$u5m.draws, by = c('region', 'strata', 'period')) %>%
  # remove strata that do not exist
  filter(!(region %in% c('mombasa', 'nairobi') & strata == 'rural')) %>% 
  # label appropriately, deselect unneeded columns, organise 
  mutate(prop = prop.x) %>%
  select(-prop.x, -prop.y) %>%
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
  select(-starts_with('theta:'))  %>%
  mutate(width = upper - lower) %>%
  select(region, period, median, lower, upper, width) %>%
  as.data.frame()

# ap
apSubnatFinal <- 
  # join the complete weights with u5m stratified results
  left_join(subnationalWeights3, apRes$u5m.draws, by = c('region', 'strata', 'period')) %>%
  # remove strata that do not exist
  filter(!(region %in% c('mombasa', 'nairobi') & strata == 'rural')) %>% 
  # label appropriately, deselect unneeded columns, organise 
  mutate(prop = prop.x) %>%
  select(-prop.x, -prop.y) %>%
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
  select(-starts_with('theta:'))  %>%
  mutate(width = upper - lower) %>%
  select(region, period, median, lower, upper, width) %>%
  as.data.frame()

# ac
acSubnatFinal <- 
  # join the complete weights with u5m stratified results
  left_join(subnationalWeights3, acRes$u5m.draws, by = c('region', 'strata', 'period')) %>%
  # remove strata that do not exist
  filter(!(region %in% c('mombasa', 'nairobi') & strata == 'rural')) %>% 
  # label appropriately, deselect unneeded columns, organise 
  mutate(prop = prop.x) %>%
  select(-prop.x, -prop.y) %>%
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
  select(-starts_with('theta:'))  %>%
  mutate(width = upper - lower) %>%
  select(region, period, median, lower, upper, width) %>%
  as.data.frame()

## map plots ----

# apc
## median
apcMedianPlot <-
  my.map.plot(data = apcSubnatFinal %>% mutate(region = region %>% as.character) %>% as.data.frame,
              geo = kenyaPoly,
              variables = 'period',
              value = 'median',
              by.data = 'region',
              by.geo = 'REGNAME',
              ncol = 3,
              per1000 = TRUE,
              direction = - 1,
              color.option = 'D',
              legend.label = 'U5MR (deaths per 1000 live births)') +
  my.map.theme(text = element_text(size = 15),
               legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
# width
apcWidthPlot <-
  my.map.plot(data = apcSubnatFinal %>% mutate(region = region %>% as.character) %>% as.data.frame,
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
  my.map.theme(text = element_text(size = 15),
               legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# ap
## median
apMedianPlot <-
  my.map.plot(data = apSubnatFinal %>% mutate(region = region %>% as.character) %>% as.data.frame,
              geo = kenyaPoly,
              variables = 'period',
              value = 'median',
              by.data = 'region',
              by.geo = 'REGNAME',
              ncol = 3,
              per1000 = TRUE,
              direction = - 1,
              color.option = 'D',
              legend.label = 'U5MR (deaths per 1000 live births)') +
  my.map.theme(text = element_text(size = 15),
               legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
# width
apWidthPlot <-
  my.map.plot(data = apSubnatFinal %>% mutate(region = region %>% as.character) %>% as.data.frame,
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
  my.map.theme(text = element_text(size = 15),
               legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# ac
## median
acMedianPlot <-
  my.map.plot(data = acSubnatFinal %>% mutate(region = region %>% as.character) %>% as.data.frame,
              geo = kenyaPoly,
              variables = 'period',
              value = 'median',
              by.data = 'region',
              by.geo = 'REGNAME',
              ncol = 3,
              per1000 = TRUE,
              direction = - 1,
              color.option = 'D',
              legend.label = 'U5MR (deaths per 1000 live births)') +
  my.map.theme(text = element_text(size = 15),
               legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
# width
acWidthPlot <-
  my.map.plot(data = acSubnatFinal %>% mutate(region = region %>% as.character) %>% as.data.frame,
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
  my.map.theme(text = element_text(size = 15),
               legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# save
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational Maps'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational Maps'))
}

ggsave(plot = apcMedianPlot,
       filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational Maps/apcMedianPlotType', type.st, '.png'),
       height = 10, width = 10)
ggsave(plot = apcWidthPlot,
       filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational Maps/apcWidthPlotType', type.st, '.png'),
       height = 10, width = 10)

ggsave(plot = apMedianPlot,
       filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational Maps/apMedianPlotType', type.st, '.png'),
       height = 10, width = 10)
ggsave(plot = apWidthPlot,
       filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational Maps/apWidthPlotType', type.st, '.png'),
       height = 10, width = 10)

ggsave(plot = acMedianPlot,
       filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational Maps/acMedianPlotType', type.st, '.png'),
       height = 10, width = 10)
ggsave(plot = acWidthPlot,
       filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational Maps/acWidthPlotType', type.st, '.png'),
       height = 10, width = 10)

## APC Estimates vs Other Estimates (subnational) ----

# data organisation
dirTemp <-
  dir %>% mutate(period = years)

subNatRes <- 
  rbind(apcSubnatFinal %>% mutate(type = haven::as_factor('APC')), 
        apSubnatFinal %>% mutate(type = haven::as_factor('AP')), 
        acSubnatFinal %>% mutate(type = haven::as_factor('AC'))) %>%
  left_join(., dirTemp, by = c('region', 'period')) %>%
  mutate(modEst = SUMMER::logit(median),
         dirEst = SUMMER::logit(mean)) %>%
  select(region, period, modEst, dirEst, type)

## direct
apcVsdirSubNat<-
  ggplot(subNatRes %>% filter(type == 'APC'), aes(x = dirEst, y = modEst)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlim(-7, -1) + ylim(-7, -1) +
  labs(x = 'Direct Estimate (logit)', y = 'APC Estimate (logit)') +
  my.theme(text = element_text(size = 15))

# plot
## ap
apcVSapSubNat<-
  ggplot() +
  geom_point(aes(y = SUMMER::logit(apcSubnatFinal$median), x = SUMMER::logit(apSubnatFinal$median))) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(-4, -2) + ylim(-4, -2) +
  labs(y = 'APC Estimate (logit)', x = 'AP Estimate (logit)') +
  my.theme(text = element_text(size = 15))

## ac
apcVSacSubNat<-
  ggplot() +
  geom_point(aes(y = SUMMER::logit(apcSubnatFinal$median), x = SUMMER::logit(acSubnatFinal$median))) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(-4, -2) + ylim(-4, -2) +
  labs(y = 'APC Estimate (logit)', x = 'AC Estimate (logit)') +
  my.theme(text = element_text(size = 15))

# save
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational APC Vs Other Estimates'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational APC Vs Other Estimates'))
}

# saving
height = width = 7

# direct
ggsave(plot = apcVsdirSubNat, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational APC Vs Other Estimates/apcVsdirSubNatType', type.st, '.png'),
       height = height, width = width)

## ap
ggsave(plot = apcVSapSubNat, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational APC Vs Other Estimates/apcVSapSubNatType', type.st, '.png'),
       height = height, width = width)
## ac
ggsave(plot = apcVSacSubNat, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Subnational APC Vs Other Estimates/apcVSacSubNatType', type.st, '.png'),
       height = height, width = width)

## spaghetti plots ----

apcSpaghettiPlot <-
  ggplot(apcSubnatFinal %>% 
           mutate(region = region %>% str_to_title,
                  median = median * 1000), 
         aes(x = period, y = median, group = region, colour = region)) +
  geom_line() +
  labs(y = 'U5MR (deaths per 1000 live births)', x= 'Period') +
  my.theme(axis.text.x = element_text(angle = 45, hjust = 1),
           strip.text.y = element_blank(),
           legend.position = 'none',
           text = element_text(size = 15))

apSpaghettiPlot <-
  ggplot(apSubnatFinal %>% 
           mutate(region = region %>% str_to_title,
                  median = median * 1000), 
         aes(x = period, y = median, group = region, colour = region)) +
  geom_line() +
  labs(y = 'U5MR (deaths per 1000 live births)', x= 'Period') +
  my.theme(axis.text.x = element_text(angle = 45, hjust = 1),
           strip.text.y = element_blank(),
           legend.position = 'none',
           text = element_text(size = 15))

acSpaghettiPlot <-
  ggplot(acSubnatFinal %>% 
           mutate(region = region %>% str_to_title,
                  median = median * 1000), 
         aes(x = period, y = median, group = region, colour = region)) +
  geom_line() +
  labs(y = 'U5MR (deaths per 1000 live births)', x= 'Period') +
  my.theme(axis.text.x = element_text(angle = 45, hjust = 1),
           strip.text.y = element_blank(),
           legend.position = 'none',
           text = element_text(size = 15))

# save
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Spaghetti Plots'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Spaghetti Plots'))
}

ggsave(apcSpaghettiPlot, file = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Spaghetti Plots/apcSpaghettiPlotType', type.st, '.png'),
       heigh = 7, width = 7)

ggsave(apSpaghettiPlot, file = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Spaghetti Plots/apSpaghettiPlotType', type.st, '.png'),
       heigh = 7, width = 7)

ggsave(acSpaghettiPlot, file = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Spaghetti Plots/acSpaghettiPlotType', type.st, '.png'),
       heigh = 7, width = 7)


# RW2 Curves ----

# data
## apc
apcCurvs <-
  rbind(apcFit$fit$summary.random$ageC_id %>% mutate(Effect = haven::as_factor('Age')),
        apcFit$fit$summary.random$perC_id %>% mutate(Effect = haven::as_factor('Period')),
        apcFit$fit$summary.random$cohC_id %>% mutate(Effect = haven::as_factor('Cohort')))
## ap
apCurvs <-
  rbind(apFit$fit$summary.random$ageC_id %>% mutate(Effect = haven::as_factor('Age')),
        apFit$fit$summary.random$perC_id %>% mutate(Effect = haven::as_factor('Period')))
# ac
acCurvs <-
  rbind(acFit$fit$summary.random$ageC_id %>% mutate(Effect = haven::as_factor('Age')),
        acFit$fit$summary.random$cohC_id %>% mutate(Effect = haven::as_factor('Cohort')))

# plots
## apc
### age
apcAgeCurvPlot<-
  ggplot(apcCurvs %>% filter(Effect == 'Age'), aes(x = ID, y = mean)) +
  geom_line(linetype = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha=0.2, linetype = 2, color = 'black') +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(y = 'Curvature', x = 'Age') +
  my.theme()
### period
apcPerCurvPlot<-
  ggplot(apcCurvs %>% filter(Effect == 'Period'), aes(x = ID, y = mean)) +
  geom_line(linetype = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha=0.2, linetype = 2, color = 'black') +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(y = 'Curvature', x = 'Period') +
  my.theme()
### curvature
apcCohCurvPlot<-
  ggplot(apcCurvs %>% filter(Effect == 'Cohort'), aes(x = ID, y = mean)) +
  geom_line(linetype = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha=0.2, linetype = 2, color = 'black') +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(y = 'Curvature', x = 'Cohort') +
  my.theme()

## ap
### age
apAgeCurvPlot<-
  ggplot(apCurvs %>% filter(Effect == 'Age'), aes(x = ID, y = mean)) +
  geom_line(linetype = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha=0.2, linetype = 2, color = 'black') +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(y = 'Curvature', x = 'Age') +
  my.theme()
### period
apPerCurvPlot<-
  ggplot(apCurvs %>% filter(Effect == 'Period'), aes(x = ID, y = mean)) +
  geom_line(linetype = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha=0.2, linetype = 2, color = 'black') +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(y = 'Curvature', x = 'Period') +
  my.theme()

## ac
### age
acAgeCurvPlot<-
  ggplot(acCurvs %>% filter(Effect == 'Age'), aes(x = ID, y = mean)) +
  geom_line(linetype = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha=0.2, linetype = 2, color = 'black') +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(y = 'Curvature', x = 'Age') +
  my.theme()
### cohort
acCohCurvPlot<-
  ggplot(acCurvs %>% filter(Effect == 'Cohort'), aes(x = ID, y = mean)) +
  geom_line(linetype = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha=0.2, linetype = 2, color = 'black') +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(y = 'Curvature', x = 'Cohort') +
  my.theme()


# save
if(!dir.exists(paths = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Temporal Curvatures'))) {
  dir.create(path = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Temporal Curvatures'))
}

# saving
height = width = 7

## apc
ggsave(plot = apcAgeCurvPlot, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Temporal Curvatures/apcAgeCurvPlotType', type.st, '.png'),
       height = height, width = width)
ggsave(plot = apcPerCurvPlot, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Temporal Curvatures/apcPerCurvPlotType', type.st, '.png'),
       height = height, width = width)
ggsave(plot = apcCohCurvPlot, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Temporal Curvatures/apcCohCurvPlotType', type.st, '.png'),
       height = height, width = width)

## ap
ggsave(plot = apAgeCurvPlot, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Temporal Curvatures/apAgeCurvPlotType', type.st, '.png'),
       height = height, width = width)
ggsave(plot = apPerCurvPlot, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Temporal Curvatures/apPerCurvPlotType', type.st, '.png'),
       height = height, width = width)

## ac
ggsave(plot = acAgeCurvPlot, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Temporal Curvatures/acAgeCurvPlotType', type.st, '.png'),
       height = height, width = width)
ggsave(plot = acCohCurvPlot, filename = paste0(kenyaPath, 'Code/Results/Plots/Type ', type.st, '/Temporal Curvatures/acCohCurvPlotType', type.st, '.png'),
       height = height, width = width)

  


