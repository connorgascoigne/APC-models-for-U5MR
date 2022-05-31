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

# data ----

# when not using region X urbanicity strata
datTemp <- getAPCdata(data = kenyaBirths, year.cut = seq(2006, 2015, by = 1), strata = c('v025'))
dat <- region.change(data = datTemp, geoPoints = kenyaPoints)

# number of observations per year ----

## monthly risk only ----

data2 <-
  dat %>%
  group_by(period) %>%
  summarise(y = sum(y),
            N = sum(N))

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


p1 <-
  ggplot2::ggplot(data2, aes(x = period, y = N)) +
  ggplot2::geom_bar(stat = 'identity') + 
  ggplot2::labs(x = 'Period', y = 'Number of Hazards') +
  # ggplot2::scale_y_continuous(label = scientific_10) +
  ggplot2::geom_text(aes(label = N), vjust = -.5) +
  my.theme(); p1

ggsave(p1, file = paste0(kenyaPath, 'Code/Results/Plots/monthly_hazards.png'),
      height = 7, width = 7)

## both month risk and observations ----

data2 <-
  dat %>%
  group_by(period) %>%
  summarise(y = sum(y),
            N = sum(N)) %>%
  pivot_longer(cols = c('y', 'N'),
               names_to = 'Type')


ggplot(data2, aes(x = period, y = value)) +
  ggplot2::geom_bar(stat = 'identity') + 
  ggplot2::facet_wrap(~Type, scale = 'free', ncol = 1) +
  ggplot2::geom_text(aes(label = value), vjust = -.5)