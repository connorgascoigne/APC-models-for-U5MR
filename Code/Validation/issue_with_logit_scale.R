# some old stuff ----
dirTrue <-
  dir %>%
  filter(years == last_year,
         region != 'All') %>%
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(years),
         true = mean) %>%
  select(region, period, true)


apcPE <-
  apcResTemp %>%
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period)) %>%
  select(region, period, median) %>%
  left_join(., dirTrue, by = c('region', 'period')) %>% 
  # relative bias for each posterior value
  mutate(bias = (median - true),
         RMSE = (median - true)^2,
         logit.bias = (SUMMER::logit(median) - SUMMER::logit(true)),
         # logit.bias = log( median*(1-true) / true*(1 - median)),
         # logit.bias = abs(SUMMER::logit(median)) - abs(SUMMER::logit(true)),
         logit.RMSE = (SUMMER::logit(median) - SUMMER::logit(true))^2,
         model = 'APC') %>%
  drop_na() %>%
  group_by(model) %>%
  summarise(bias = mean(bias),
            RMSE = sqrt(mean(RMSE)),
            logit.bias = mean(logit.bias),
            logit.RMSE = sqrt(mean(logit.RMSE)))

apPE <-
  apResTemp %>%
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period)) %>%
  select(region, period, median) %>%
  left_join(., dirTrue, by = c('region', 'period')) %>%
  # relative bias for each posterior value
  mutate(bias = (median - true),
         RMSE = (median - true)^2,
         logit.bias = (SUMMER::logit(median) - SUMMER::logit(true)),
         # logit.bias = log( median*(1-true) / true*(1 - median)),
         # logit.bias = abs(SUMMER::logit(median)) - abs(SUMMER::logit(true)),
         logit.RMSE = (SUMMER::logit(median) - SUMMER::logit(true))^2,
         model = 'AP') %>%
  drop_na() %>%
  group_by(model) %>%
  summarise(bias = mean(bias),
            RMSE = sqrt(mean(RMSE)),
            logit.bias = mean(logit.bias),
            logit.RMSE = sqrt(mean(logit.RMSE)))

# logit(apPE$median) - logit(apPE$true) == apPE$logit.bias
# (logit(apPE$median) - logit(apPE$true))^2 == apPE$logit.RMSE

acPE <-
  acResTemp %>%
  mutate(region = haven::as_factor(region),
         period = haven::as_factor(period)) %>%
  select(region, period, median) %>%
  left_join(., dirTrue, by = c('region', 'period')) %>% 
  # relative bias for each posterior value
  mutate(bias = (median - true),
         RMSE = (median - true)^2,
         logit.bias = (SUMMER::logit(median) - SUMMER::logit(true)),
         # logit.bias = log( median*(1-true) / true*(1 - median)),
         # logit.bias = abs(SUMMER::logit(median)) - abs(SUMMER::logit(true)),
         logit.RMSE = (SUMMER::logit(median) - SUMMER::logit(true))^2,
         model = 'AC') %>%
  drop_na() %>%
  group_by(model) %>%
  summarise(bias = mean(bias),
            RMSE = sqrt(mean(RMSE)),
            logit.bias = mean(logit.bias),
            logit.RMSE = sqrt(mean(logit.RMSE)))

peScores <- 
  rbind(apcPE, apPE, acPE) %>%
  mutate(bias = bias * 1000,
         RMSE = RMSE * 100,
         logit.bias = logit.bias * 1000,
         logit.RMSE = logit.RMSE * 100)

print(xtable::xtable(peScores, type = 'latex'), include.rownames = FALSE, digits = 2)
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


