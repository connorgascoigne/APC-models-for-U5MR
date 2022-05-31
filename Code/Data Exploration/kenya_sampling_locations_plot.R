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

# plotting clusters ----


# organise spatial locations
geo1 <- 
  ggplot2::fortify(kenyaPoly, region = 'REGNAME') %>%
  arrange(order)

urbanCluster <-
  kenyaBirths %>%
  select(v001, v025) %>%
  rename(cluster = v001, urban = v025) %>%
  mutate(urban = factor(urban, labels = c('Urban', 'Rural'))) %>%
  unique()

clustPoints <-
  kenyaPoints %>%
  as.data.frame() %>%
  mutate(id = ADM1NAME %>% tolower()) %>%
  select(DHSCLUST, id, LONGNUM, LATNUM) %>%
  filter(LONGNUM > 30) %>%
  rename(cluster = DHSCLUST, long = LONGNUM, lat = LATNUM) %>%
  left_join(., urbanCluster, by = 'cluster') %>%
  na.omit

# colour
p <-
  ggplot2::ggplot(clustPoints, aes(x = long, y = lat, color = urban)) +
  ggplot2::geom_polygon(data = geo1, aes(x = long, y = lat, group = group), fill = 'white', color = 'gray20', size = .5) +
  ggplot2::coord_map() +
  ggplot2::geom_point(aes(shape = urban), size = 1) +
  ggplot2::scale_shape_manual(values=c(8, 0)) +
  ggplot2::labs(x = 'Long', y = 'Lat') +
  my.map.theme(text = element_text(size = 15),
               legend.position = 'top',
               legend.title = element_blank()); p

# black and white
# p <- 
#   ggplot2::ggplot(clustPoints, aes(x = long, y = lat)) +
#   ggplot2::geom_polygon(data = geo1, aes(x = long, y = lat, group = group), fill = 'white', color = 'gray20', size = .5) + 
#   ggplot2::coord_map() +
#   ggplot2::geom_point(aes(shape = urban), size = 1) + 
#   ggplot2::scale_shape_manual(values=c(8, 0)) +
#   ggplot2::labs(x = 'Long', y = 'Lat') +
#   my.map.theme(legend.position = 'top',
#                legend.title = element_blank()); p

ggsave(file = paste0(kenyaPath, 'Code/Plots/kenya_sample_locations.png'),
       plot = p,
       height = 7, width = 7)

# ggsave(file = 'kenya_sample_locations.png',
#        plot = p,
#        height = 7, width = 7)


