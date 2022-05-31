# packages ----

library(haven) # data import
library(tidyverse) # just bare stuff for R
library(SUMMER) # U5MR package from wakeyyy
library(foreign) # importing data
library(INLA) # fitting models
library(rgdal) # spatial modelling
library(patchwork)

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


## model fit ----

apFit <-
  reparam.model(data = dat, mod = 'ap',
                inla.mode = 'experimental', is.strata = TRUE, intercept = FALSE, Amat = NULL, 
                control.compute = list(config = TRUE), safe = TRUE)

## posterior summary ----

apRes <- getSamples(apFit, strata.weights = nationalWeights, save.u5m.draws = TRUE)

## smooth cluster ----

smthClusRes <- readRDS(file = paste0(kenyaPath, 'Code/Results/SUMMER/smooth_cluster_res.rds'))

## subnational AP ----

apSubRes <- readRDS(file = paste0(kenyaPath, 'Code/Results/Type 4/ap/model_res.rds'))

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

apSubFinal <- 
  # join the complete weights with u5m stratified results
  left_join(subnationalWeights2, apSubRes$u5m.draws, by = c('region', 'strata', 'period')) %>%
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
  mutate(type = haven::as_factor('AP Subnational')) %>%
  select(region, period, median, lower, upper, type) %>%
  as.data.frame()

# plot ----

## without subnational ----

natRes <-
  rbind(smthClusRes$overall %>%
          mutate(period = years,
                 type = haven::as_factor('AT National')) %>%
          select(region, period, median, lower, upper, type),
        apRes$overall %>%
          mutate(type = haven::as_factor('AP National')) %>%
          select(region, period, median, lower, upper, type))

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

natPlotError + natPlotLine

ggsave(natPlotError, file = paste0(kenyaPath, 'Code/Results/Plots/Age Time SUMMER/natPlotErrorSUMMER.png'), height = 7, width = 7)
ggsave(natPlotLine, file = paste0(kenyaPath, 'Code/Results/Plots/Age Time SUMMER/natPlotLineSUMMER.png'), height = 7, width = 7)

## with subnational ----

natRes <-
  rbind(smthClusRes$overall %>%
          mutate(period = years,
                 type = haven::as_factor('AT National')) %>%
          select(region, period, median, lower, upper, type),
        apRes$overall %>%
          mutate(type = haven::as_factor('AP National')) %>%
          select(region, period, median, lower, upper, type),
        apSubFinal %>%
          mutate(type = haven::as_factor('AP Subnational')))

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

natPlotError + natPlotLine

ggsave(natPlotError, file = paste0(kenyaPath, 'Code/Results/Plots/Age Time SUMMER/natPlotErrorSUMMER (with subnational).png'), height = 7, width = 7)
ggsave(natPlotLine, file = paste0(kenyaPath, 'Code/Results/Plots/Age Time SUMMER/natPlotLineSUMMER (with subnational).png'), height = 7, width = 7)
