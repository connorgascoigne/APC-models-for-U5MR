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

# plots of sampling proportions ---- 

# organise spatial locations
geo1 <- 
  ggplot2::fortify(kenyaPoly, region = 'REGNAME') %>%
  arrange(order)

# combine geo data and model results
data2 <- 
  merge(geo1, subnationalWeights, by = 'id', by.y = 'region') %>% # need sort = FALSE here to stop polygon ordering being messed up
  arrange(order) # arrange by the fortifies objects order

strataPropPlot <- 
  ggplot2::ggplot(data2 %>% filter(period %in% c(2006, 2010, 2014)) %>% mutate(urban = str_to_title(urban))) +
  ggplot2::geom_polygon(aes(x = long, y = lat, group = group, fill = prop), color = 'gray20', size = .5) +
  ggplot2::scale_fill_viridis_c('Strata Proportions', direction = -1, option = 'D') +
  ggplot2::coord_map() +
  ggplot2::facet_grid(urban~period, switch = 'y') + 
  labs(x = 'Long', y = 'Lat') +
  my.map.theme(legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

strataPropPlot

data3 <-
  merge(geo1, admin1Proportions, by = 'id', by.y = 'region')

regionPropPlot <-
  ggplot2::ggplot(data3 %>% filter(period %in% c(2006, 2010, 2014))) +
  ggplot2::geom_polygon(aes(x = long, y = lat, group = group, fill = prop), color = 'gray20', size = .5) +
  ggplot2::scale_fill_viridis_c('Region Proportions', direction = -1, option = 'A') +
  ggplot2::coord_map() +
  ggplot2::facet_wrap(~period, ncol = 3) + 
  labs(x = 'Long', y = 'Lat') +
  my.map.theme(legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

regionPropPlot
  
ggsave(file = paste0(kenyaPath, 'Code/Plots/strata_proportions_plot.png'),
       plot = strataPropPlot,
       height = 7, width = 7)

ggsave(file = paste0(kenyaPath, 'Code/Plots/region_proportions_plot.png'),
       plot = regionPropPlot,
       height = 7, width = 7)

