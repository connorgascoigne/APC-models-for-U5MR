# organise data for APC
getAPCdata <- function(data, year.cut = seq(1980, 2020, by = 1), age.cut = c(0, 1, 12, 24, 36, 48, 60), strata = c('v024', 'v025')){
  
  # data = kenyaBirths; year.cut = seq(1980, 2020, by = 5); age.cut = c(0, 1, 12, 24, 36, 48, 60); strata = c('v024', 'v025')
  
  # get birth data
  datTemp <- SUMMER::getBirths(data = data, month.cut = age.cut, year.cut = year.cut, strata = strata)
  
  # Suppress summarise info
  options(dplyr.summarise.inform = FALSE)
  
  # organise temporal
  ## columns include age, time (period), dob (cohort)
  datTemp2 <- 
    datTemp[,c('v001', 'v005', 'v024', 'v025', 'strata', 'age', 'time', 'dob', 'died')] %>%
    mutate(dob = floor((dob-1)/12)+1900) %>%
    group_by(v001, v005, v024, v025, strata, age, time, dob) %>%
    summarise(died = sum(died),
              total = n()) %>%
    as.data.frame()
  
  # add column for age bin and age bin midpoint
  ## setting up
  if (age.cut[1] == 0) 
    age.cut <- age.cut[-1]
  ageBinsMid <- rep(NA, length(age.cut))
  ageBinsMid[1] <- 0.5
  ## make unique vector of bin and bin midpoint
  if (length(age.cut) > 1) {
    for (i in 1:(length(age.cut) - 1)) {
      ageBinsMid[i + 1] <- mean(c(age.cut[i + 1] - 1, age.cut[i]))
    }
  }
  ## go down data frame and add additional columns for bin and bin midpoint
  datTemp2$ageMid <- ageBinsMid[1]
  ageBins <- datTemp2$age %>% unique %>% sort %>% as.vector
  for (i in 1:(length(age.cut) - 1)) {
    datTemp2$ageMid[datTemp2$age == ageBins[i + 1]] <- ageBinsMid[i + 1]
  }
  
  perBinsMid <- rep(NA, (length(year.cut)-1))
  ## make unique vector of bin and bin midpoint
  for (i in 1:(length(year.cut) - 1)) {
    perBinsMid[i] <- mean(c(year.cut[i + 1], year.cut[i]))
  }
  ## go down data frame and add additional columns for bin and bin midpoint
  datTemp2$perMid <- perBinsMid[1]
  perBins <- datTemp2$time %>% unique %>% sort %>% as.vector
  for (i in 1:(length(year.cut) - 1)) {
    datTemp2$perMid[datTemp2$time == perBins[i + 1]] <- perBinsMid[i + 1]
  }
  
  datTemp3 <-
    datTemp2[,c('v001', 'v005', 'v024', 'v025', 'strata', 'age', 'ageMid', 'time', 'perMid', 'dob', 'died', 'total')] %>%
    mutate(v005 = v005/1e6, # DHS recode manual
           v024 = haven::as_factor(v024),
           v025 = haven::as_factor(v025),
           strata = haven::as_factor(strata)) %>%
    rename(cluster = v001, weights = v005, region = v024, urban = v025, period = time, cohort = dob, y = died, N = total) %>%
    group_by(cluster, weights, region, urban, age, period, cohort) %>%
    arrange %>%
    as.data.frame
  
  datTemp3
  
}

# update the regions
region.change <- function(data, geoPoints){
  
  spObj <- data.frame(cluster = geoPoints$DHSCLUST, region = geoPoints$ADM1NAME, lon = geoPoints$LONGNUM, lat = geoPoints$LATNUM)
  idx <- match(data$cluster, spObj$cluster)
  dataNew <-
    data %>%
    mutate(region = spObj$region[idx] %>% tolower %>% haven::as_factor())
  
  dataNew
  
}

# reparameterise model
reparam.model <- function(data, mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')[1], Amat = NULL,
                          slopeDrop = NULL, family = c('betabinomial', 'binomial')[1],
                          is.strata = FALSE, intercept = FALSE,
                          hyper = c('pc', 'gamma')[1], formula = NULL,
                          pc.u = 1, pc.alpha = 0.01, pc.u.phi = 0.5, pc.alpha.phi = 2/3,
                          overdisp.mean = 0, overdisp.prec = 0.4, 
                          type.st = 4, pc.st.u = NA, pc.st.alpha = NA,
                          control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                          inla.mode = c('classic', 'twostage', 'experimental')[1],
                          control.compute = list(config = TRUE), verbose = FALSE, ...){
  
  # data = dat; mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')[1]; Amat = Amat;
  # slopeDrop = 'c'; family = c('betabinomial', 'binomial')[2]
  # is.strata = TRUE
  # intercept = FALSE
  # hyper = c('pc', 'gamma')[1]
  # formula = NULL
  # pc.u = 1; pc.alpha = 0.01; pc.u.phi = 0.5; pc.alpha.phi = 2/3; overdisp.mean = 0; overdisp.prec = 0.4;
  # type.st = 1; pc.st.u = NA; pc.st.alpha = NA;
  # control.inla = list(strategy = 'adaptive', int.strategy = 'auto'); control.compute = list(config = TRUE); inla.mode = c('classic', 'twostage', 'experimental')[3]; verbose = FALSE
  
  if (!"config" %in% names(control.compute)) {
    message("config = TRUE is added to control.compute so that posterior draws can be taken.")
    control.compute$config <- TRUE
  }
  if(mod == 'apc' && is.null(slopeDrop)) stop('Warning: When fitting an APC model need to drop linear column of one of age, period or cohort.')
  if(mod == 'apc' &&! (slopeDrop %in% c('a', 'p', 'c'))) stop('slopeDrop needs to be one of a, p or c')
  
  # will this be a spatial model or not
  if (!is.null(Amat)) {
    if (is.null(rownames(Amat))) {
      stop('Row names of Amat needs to be specified to region names.')
    }
    if (is.null(colnames(Amat))) {
      stop('Column names of Amat needs to be specified to region names.')
    }
    if (sum(rownames(Amat) != colnames(Amat)) > 0) {
      stop('Row and column names of Amat needs to be the same.')
    }
    is.spatial <- TRUE
  } else {
    # Amat <- matrix(1, 1, 1)
    # colnames(Amat) <- rownames(Amat) <- 'All'
    data$region <- "All" %>% as.factor
    is.spatial <- FALSE
  }
  
  # data augmentation
  ## scale the temporal terms
  ## define the additional columns for curvature
  ## make numeric labels for region, urban and strata and scale
  ### scale centers and then scales between 0 and 1 
  dataAug <- 
    data  %>%
    # complete so all region/strata combiantions have each APC combination
    tidyr::complete(nesting(region, urban, strata),
                    nesting(age, ageMid, period, perMid, cohort)) %>%
    # labels for INLA modelling
    mutate(space_time = interaction(region, period),
           cluster_id = cluster %>% as.factor %>% as.numeric,
           region_id = region %>% as.numeric,
           region_struct_id = region_id,
           space_time_id = space_time %>% as.numeric,
           urban_id = urban %>% as.numeric,
           strata_id = strata %>% as.numeric,
           age_id = ageMid,
           per_id = perMid,
           coh_id = cohort,
           ageC_id = age_id, 
           perC_id = per_id, 
           cohC_id = coh_id)
  
  # define the formula for INLA
  if (is.null(formula)) {
    
    if(is.spatial){
      
      # number of periods and regions
      P <- dataAug$per_id %>% unique %>% length
      R <- dataAug$region_id %>% unique %>% length
      
      # precision matrices
      ## time
      ### unstructured
      Q_time_unstruc <- diag(P)
      ### structured
      Q_time_struc <- INLA:::inla.rw(n = P, order = 2, scale.model = TRUE, sparse = TRUE)
      ## space
      ### unstructured
      Q_space_unstruc <- diag(R)
      ### structured
      Q_space_struc <- Amat
      #### organise Amat into ICAR precision
      if (sum(Q_space_struc > 0 & Q_space_struc < 1) != 0) {
        #### if 0 < Q_{ij} < 1 then make Q_{ij} = 1
        for (i in 1:nrow(Q_space_struc)) {
          idx <- which(Q_space_struc[i, ] > 0 & Q_space_struc[i, ] < 1)
          Q_space_struc[i, idx] <- 1
        }
      }
      #### turn Amat into an ICAR precision
      ##### no 0s on diagonal; main diag to be the number of neighbours in total; off diagonal to be 0 for not neighbours and -1 for neighbours
      diag(Q_space_struc) <- 0
      diag <- apply(Q_space_struc, 1, sum)
      Q_space_struc[Q_space_struc != 0] <- -1
      diag(Q_space_struc) <- diag
      Q_space_struc <- INLA:::inla.scale.model(Q_space_struc, list(A = matrix(1, 1, dim(Q_space_struc)[1]), e = 0))
    }
    
    if (tolower(hyper) == 'pc') {
      
      # pc hyper priors
      ## temproal
      hyper_pc_time <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
      ## spatial
      hyper_pc_space_struc <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)), phi = list(prior = 'pc', param = c(pc.u.phi, pc.alpha.phi)))
      # hyper_pc_space_struc <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)), theta = list(prior = 'pc', param = c(pc.u.phi, pc.alpha.phi)))
      hyper_pc_space_unstruc <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
      ## spatio-temproal
      pc.st.u <- ifelse(is.na(pc.st.u), pc.u, pc.st.u)
      pc.st.alpha <- ifelse(is.na(pc.st.alpha), pc.alpha, pc.st.alpha)
      hyper_pc_space_time <- list(prec = list(prior = "pc.prec", param = c(pc.st.u, pc.st.alpha)))
      
      # full temporal formula
      formula <- 
        y ~ age_id + per_id + coh_id +
        f(ageC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2) +
        f(perC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2) +
        f(cohC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2)
      
      message("----------------------------------",
              "\nModel Specification", 
              "\n Temporal effect(s):               ",
              "\n  APC type:                    ",
              mod,
              appendLF = FALSE)
      # update pre and post forumla for the temporal models
      # # do not remove any linear terms from pre fit as they are all needed later in 
      # # re-parameterisation
      if (mod == 'apc'){
        message("\n  Slope dropped:                 ", slopeDrop, appendLF = FALSE)
        if(slopeDrop == 'a'){
          formula <- update(formula, ~. - age_id)
        } else if (slopeDrop == 'p'){
          formula <- update(formula, ~. - per_id)
        } else if (slopeDrop == 'c'){
          formula <- update(formula, ~. - coh_id)
        }
      } else if (mod == 'ac'){
        formula <- update(formula, ~.
                          - per_id - f(perC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2))
      } else if (mod == 'pc'){
        formula <- update(formula, ~.
                          - age_id - f(ageC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2))
      } else if (mod == 'ap'){
        formula <- update(formula, ~.
                          - coh_id - f(cohC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2))
      } else if (mod == 'c'){
        formula <- update(formula, ~.
                          - age_id - f(ageC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2) -
                            per_id - f(perC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2))
      } else if (mod == 'a'){
        formula <- update(formula, ~.
                          - per_id - f(perC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2) -
                            coh_id - f(cohC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2))
      } else if(mod == 'p'){
        formula <- update(formula, ~.
                          - age_id - f(ageC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2) -
                            coh_id - f(cohC_id, model = 'rw2', hyper = hyper_pc_time, constr = TRUE, rankdef = 2))
      }
      
      if (is.spatial) {
        message("\n Spatial effect(s):               ",
                "\n  Structured:                 bym2",
                "\n  Interaction type:              ", type.st,
                appendLF = FALSE)
        formula <- 
          update(formula, ~. +
                   f(region_struct_id, graph = Amat, model = 'bym2', hyper= hyper_pc_space_struc, scale.model = TRUE, adjust.for.con.comp = TRUE))
        
        
        # space-time interaction
        ## defining precision and constraint matrices
        ## order matters here! In our data set, the time ordering take precedent over the space ordering
        ## so we must have the time Q on the left and space Q on the right in kronecker products
        if (type.st == 1) {
          # kronecker product for the Type I Interaction
          ## IID x IID
          formula <- update(formula, ~. + f(space_time_id, model="iid", hyper = hyper_pc_space_time))
        } else {
          
          if (type.st == 2) {
            # kronecker product for the Type II Interaction
            ## RW2 x IID
            Q_space_time <- Q_time_struc %x% Q_space_unstruc
          } else if (type.st == 3) {
            # kronecker product for the Type III Interaction
            ## IID x ICAR
            Q_space_time <- Q_time_unstruc %x% Q_space_struc
          } else {
            # kronecker product for the Type IV Interaction
            ## RW2 x ICAR
            Q_space_time <- Q_time_struc %x% Q_space_struc
          }
          
          # constraints for space-time term
          ## sum to zero constraints for identifiability
          ## rank def is due to 0 eigenvalues
          ## define constraint matrix from the eigenvectors relating to zero eigenvalues
          ## rank def is the number of 0 eigenvalues
          eigenQ <- eigen(Q_space_time)
          ids <- which(eigenQ$values < 1e-10)
          cMat <- t(eigenQ$vectors[,ids])
          
          # # check the dimensions line up correctly
          # R*(P-2) + nrow(cMat) == nrow(Q_space_time) # Type II
          # (R-1)*P + nrow(cMat) == nrow(Q_space_time) # Type III
          # (R-1)*(P-2) + nrow(cMat) == nrow(Q_space_time) # Type IV
          
          
          # add spatio-temporal term to formula
          formula <- update(formula, ~. +
                              f(space_time_id,
                                model = "generic0",
                                Cmatrix = Q_space_time,
                                extraconstr = list(A = cMat, e = rep(0, nrow(cMat))),
                                rankdef = nrow(cMat),
                                hyper = hyper_pc_space_time))
        } 
      }
    } else if (tolower(hyper) == 'gamma') {
      
      # gamma priors
      priors <- 
        SUMMER::simhyper(R = 2, # # 95% prior interval for the residual odds ratio
                         nsamp = 1e+05, # sample to simulate for scaling factor
                         nsamp.check = 5000, # sample to simulate for checking range
                         Amat = Amat, # adjacenecy matrix of the areas in the data
                         nperiod = dataAug$period %>% unique %>% length, # numerical values of how many time periods in the data
                         only.iid = TRUE # indicator for IID or not
        )
      a.iid <- priors$a.iid
      b.iid <- priors$b.iid
      a.rw <- priors$a.iid
      b.rw <- priors$b.iid
      a.icar <- priors$a.iid
      b.icar <- priors$b.iid
      
      # full temporal formula
      formula <- 
        y ~ age_id + per_id + coh_id +
        f(ageC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2) +
        f(perC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2) +
        f(cohC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2)
      
      message("----------------------------------",
              "\nModel Specification", 
              "\n Temporal effect(s):               ",
              "\n  APC type:                    ",
              mod,
              appendLF = FALSE)
      # update pre and post forumla for the temporal models
      # # do not remove any linear terms from pre fit as they are all needed later in 
      # # re-parameterisation
      if (mod == 'apc'){
        message("\n  Slope dropped:                 ", slopeDrop, appendLF = FALSE)
        if(slopeDrop == 'a'){
          formula <- update(formula, ~. - age_id)
        } else if (slopeDrop == 'p'){
          formula <- update(formula, ~. - per_id)
        } else if (slopeDrop == 'c'){
          formula <- update(formula, ~. - coh_id)
        }
      } else if (mod == 'ac'){
        formula <- update(formula, ~.
                          - per_id - f(perC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2))
      } else if (mod == 'pc'){
        formula <- update(formula, ~.
                          - age_id - f(ageC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2))
      } else if (mod == 'ap'){
        formula <- update(formula, ~.
                          - coh_id - f(cohC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2))
      } else if (mod == 'c'){
        formula <- update(formula, ~.
                          - age_id - f(ageC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2) -
                            per_id - f(perC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2))
      } else if (mod == 'a'){
        formula <- update(formula, ~.
                          - per_id - f(perC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2) -
                            coh_id - f(cohC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2))
      } else if(mod == 'p'){
        formula <- update(formula, ~.
                          - age_id - f(ageC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2) -
                            coh_id - f(cohC_id, model = 'rw2', param = c(a.rw, b.rw), constr = TRUE, rankdef = 2))
      }
      
      if(is.spatial) {
        message("\n Spatial effect(s):               ",
                "\n  Structured:                 bym2",
                "\n  Interaction type:              ", type.st,
                appendLF = FALSE)
        formula <- 
          update(formula, ~.
                 + f(region_struct_id, graph = Amat, model = 'bym2', param = c(a.icar, b.icar), scale.model = TRUE, adjust.for.con.comp = TRUE))
        
        # space-time interaction
        ## order matters here! In our data set, the time ordering take precedent over the space ordering
        ## so we must have the time Q on the left and space Q on the right in kronecker products 
        if (type.st == 1) {
          # kronecker product for the Type I Interaction
          ## IID x IID
          formula <- update(formula, ~. + f(space_time_id, model="iid", param = c(a.iid, b.iid)))
        } else {
          if (type.st == 2) {
            # kronecker product for the Type II Interaction
            ## RW2 x IID
            Q_space_time <- Q_time_struc %x% Q_space_unstruc
          } else if (type.st == 3) {
            # kronecker product for the Type III Interaction
            ## IID x ICAR
            Q_space_time <- Q_time_unstruc %x% Q_space_struc
          } else {
            # kronecker product for the Type IV Interaction
            ## RW2 x ICAR
            Q_space_time <- Q_time_struc %x% Q_space_struc
          }
          # constraints for space-time term
          ## sum to zero constraints for identifiability
          ## rank def is due to 0 eigenvalues
          ## define constraint matrix from the eigenvectors relating to zero eigenvalues
          ## rank def is the number of 0 eigenvalues
          eigenQ <- eigen(Q_space_time)
          ids <- which(eigenQ$values < 0.00001)
          cMat <- t(eigenQ$vectors[,ids])
          
          # add spatio-temporal term to formula
          formula <- update(formula, ~. +
                              f(space_time_id,
                                model = "generic0", 
                                Cmatrix = Q_space_time, 
                                extraconstr = list(A = cMat, e = rep(0, nrow(cMat))),
                                rankdef = nrow(cMat),
                                param = c(a.iid, b.iid)))
        }
        
      }
      
    } else {
      stop('hyper needs to be either pc or gamma.')
    }
    
    message("\n Other effect(s):", appendLF = FALSE)
    
    # cluster random effect for binomial model
    if (family == 'binomial') {
      if (tolower(hyper) == 'gamma') {
        message("\n  cluster:                     iid", appendLF = FALSE)
        formula <- update(formula, ~. +  f(cluster_id, model = 'iid', param = c(a.iid, b.iid)))
      } else if (tolower(hyper) == 'pc') {
        message("\n  cluster:                     iid", appendLF = FALSE)
        formula <- update(formula, ~. +  f(cluster_id, model = 'iid', hyper = hyper_pc_time))
      } else {
        stop('hyper needs to be either pc or gamma.')
      }
    }
    
    # sorting out the fixed effects
    ## urban/rural strata
    if(is.strata){
      message("\n  Strata term included:        yes", appendLF = FALSE)
      # update the formula
      formula <- update(formula, ~. + strata_id)
    } else {
      message("\n  Strata term included:         no", appendLF = FALSE)
      dataAug$strata <- as.factor('All')
      dataAug$strata_id <- dataAug$strata_id <- 1
    }
    
    ## intercept
    if(intercept){
      message("\n  Intercept term included:     yes", appendLF = FALSE)
    } else {
      message("\n  Intercept term included:      no", appendLF = FALSE)
      formula <- update(formula, ~. - 1)
    }
    
  }
  
  message("\nFitting model...",
          appendLF = FALSE)
  startClock <- proc.time() # Start the clock!
  # the extra information depending on the model
  if(family == 'betabinomial'){
    control.family <- list(hyper = list(rho = list(param = c(overdisp.mean, overdisp.prec))))
    fit <-
      INLA::inla(formula, family = family, data = dataAug, Ntrials = dataAug$N,
                 control.compute = control.compute,
                 control.family = control.family,
                 control.predictor = list(compute = FALSE, link = 1),
                 control.inla = control.inla,
                 lincomb = NULL,
                 inla.mode = inla.mode,
                 verbose = verbose,
                 ...)
  } else {
    fit <-
      INLA::inla(formula, family = family, data = dataAug, Ntrials = dataAug$N,
                 control.compute = control.compute,
                 control.predictor = list(compute = FALSE, link = 1),
                 control.inla = control.inla,
                 lincomb = NULL,
                 inla.mode = inla.mode,
                 verbose = verbose,
                 ...)
  }
  
  endClock <- proc.time() - startClock # Stop the clock
  message("\n Time take to fit model(s):     ", endClock[3] %>% as.numeric %>% round(),
          "\n----------------------------------",
          appendLF = FALSE)
  
  priors <- list(pc.u = pc.u, pc.alpha = pc.alpha, pc.u.phi = pc.u.phi, 
                 pc.alpha.phi = pc.alpha.phi, 
                 pc.st.u = pc.st.u, 
                 pc.st.alpha = pc.st.alpha, 
                 overdisp.mean = overdisp.mean, overdisp.prec = overdisp.prec)
  
  return(list(model = formula, fit = fit, family = family, 
              Amat = Amat, newdata = dataAug, type.st = type.st,
              priors = priors))
}

# posterior sampling
getSamples <- function(inla_mod, nsim = 1000, use.approximation = TRUE, verbose = FALSE, CI = 0.95, save.u5m.draws = FALSE, strata.weights = NULL, ...){
  
  # inla_mod = apcBinFit;  strata.weights = NULL; use.approximation = TRUE
  # nsim = 10; verbose = FALSE; CI = 0.95; save.u5m.draws = TRUE
  # strata.weights = NULL
  
  # lower and upper CI calculators
  lowerCI <- (1 - CI)/2
  upperCI <- 1 - lowerCI
  
  # contents tags
  cs <- inla_mod$fit$misc$configs$contents$tag
  ## remove the tags of predictor and clustid to speed up sampling
  cs <- cs[cs != "Predictor"]
  cs <- cs[cs != "cluster_id"]
  select <- list()
  for (i in 1:length(cs)) {
    select[[i]] <- 0
    names(select)[i] <- cs[i]
  }
  
  message("-----------------------------------",
          "\nStarting posterior sampling...   ",
          appendLF = F)
  start1 <- proc.time()
  sampAll <- INLA::inla.posterior.sample(n = nsim, result = inla_mod$fit, intern = TRUE, selection = select, verbose = verbose, ...)
  end1 <- proc.time() - start1
  message("\nCleaning up results...           ",
          appendLF = F)
  
  start2 <- proc.time()
  
  # to get the names for labelling
  fields <- rownames(sampAll[[1]]$latent)
  
  # data
  data2 <-
    # data frame from INLA
    inla_mod$fit$.args$data %>%
    # select the columns from original data going into function EXCEPT cluster
    select(region, urban, strata, age, ageMid, period, perMid, cohort) %>%
    # only distinct columns of the above 
    distinct() %>%
    # recreat the id variables
    mutate(region_id = region %>% as.numeric,
           region_struct_id = region_id,
           space_time_id = interaction(region, period),
           urban_id = urban %>% as.numeric,
           strata_id = strata %>% as.numeric,
           age_id = ageMid,
           per_id = perMid,
           coh_id = cohort,
           ageC_id = age_id, 
           perC_id = per_id, 
           cohC_id = coh_id)
  
  # locations
  data2.loc <- 
    data2 %>%
    # chnage the variable names to match INLA output names
    mutate(ageC_id = age_id %>% as.factor %>% as.numeric, 
           perC_id = per_id %>% as.factor %>% as.numeric, 
           cohC_id = coh_id %>% as.factor %>% as.numeric,
           region_struct_id = region_struct_id %>% as.factor %>% as.numeric,
           space_time_id = space_time_id %>% as.numeric,
           intercept = '(Intercept):1',
           strata_id = 'strata_id:1',
           age_id = 'age_id:1',
           per_id = 'per_id:1',
           coh_id = 'coh_id:1') %>% 
    # match the names to INLA output positions
    mutate(intercept = match(intercept, fields),
           strata_id = match(strata_id, fields),
           age_id = match(age_id, fields),
           per_id = match(per_id, fields),
           coh_id = match(coh_id, fields),
           ageC_id = match(paste0('ageC_id:', ageC_id), fields),
           perC_id = match(paste0('perC_id:', perC_id), fields),
           cohC_id = match(paste0('cohC_id:', cohC_id), fields),
           region_struct_id = match(paste0('region_struct_id:', region_struct_id), fields),
           space_time_id = match(paste0('space_time_id:', space_time_id), fields))
  
  # only the random effect positions (not iid ones though)
  data2.loc.random.effects <- data2.loc[,c('ageC_id', 'perC_id', 'cohC_id', 'region_struct_id', 'space_time_id')] 
  
  ## log precision of the random effect to marginalize over
  tau <- rep(NA, nsim)
  ## linear predictor
  theta <- matrix(0, nrow = nsim, ncol = dim(data2)[1])
  
  for (i in 1:nsim) {
    
    # extract sample i 
    draw <- sampAll[[i]]$latent
    
    # linear predictor
    ## drop the columns for the fixed effects
    ## only want samples from the random effect columns
    theta[i, ] <- apply(data2.loc.random.effects, 1, function(x, ff) { sum(ff[x], na.rm = TRUE) }, draw)
    
    # now include the fixed effects
    ## want to have the slope multiplied by the value
    add.intercept <- draw[data2.loc$intercept]
    add.strata.slope <- draw[data2.loc$strata_id] * data2$strata_id
    add.age.slope <- draw[data2.loc$age_id] * data2$age_id
    add.period.slope <- draw[data2.loc$per_id] * data2$per_id
    add.cohort.slope <- draw[data2.loc$coh_id] * data2$coh_id
    
    ## if there is any NA return give 0
    ### NA will be return for the linear slope that is dropped
    add.intercept[is.na(add.intercept)] <- 0
    add.strata.slope[is.na(add.strata.slope)] <- 0
    add.age.slope[is.na(add.age.slope)] <- 0
    add.period.slope[is.na(add.period.slope)] <- 0
    add.cohort.slope[is.na(add.cohort.slope)] <- 0
    
    theta[i, ] <- theta[i, ] + add.intercept + add.strata.slope + add.age.slope + add.period.slope + add.cohort.slope
    if (inla_mod$family == 'binomial') {
      tau[i] <- exp(sampAll[[i]]$hyperpar[["Log precision for cluster_id"]]) # tau = exp(log(precision))
    }
  }
  
  # function for summary
  ## quantiles, mean and varaince
  my.summary <- function(x, lowerCI, upperCI) {
    
    # x = theta1[,1]; lowerCI = lowerCI; upperCI = upperCI
    
    qntl <- quantile(x, probs = c(lowerCI, 0.5, upperCI))
    data.frame(mean = mean(x), variance = var(x), lower = qntl[1], median = qntl[2], upper = qntl[3])
    
  }
  
  # marginalisation function
  ## for binomial distribution
  ## use approximation or numerical integration
  marginalise.logit.binomial = function(tauTheta, use.approximation = TRUE, ...) {
    
    # tauTheta = cbind(tau, theta)[1,]; use.approximation = FALSE
    
    # tau = 1/sigma^2
    sigma <- sqrt(1/tauTheta[1])
    theta <- as.matrix(tauTheta[-1], ncol = 1)
    
    # calculation marginalistion
    if(sigma == 0){
      SUMMER::expit(theta)
    } else {
      if(any(is.na(c(theta, sigma)))){
        NA
      } else if(!use.approximation) {
        # numerically calculate the mean
        ## integrand to calculate numerical marginal
        integrand <- function(x, theta, sigma){ exp(plogis(x, log.p=TRUE) + dnorm(x, mean = theta, sd = sigma, log=TRUE)) }
        ## function to find integral
        ### requires theta and sigma to be single elements
        my.integrate <- function(fun, theta, sigma, ...){
          integrate(fun, lower = theta-10*sigma, upper = theta+10*sigma, theta = theta, sigma = sigma, abs.tol = 0, ...)$value
        }
        ## apply on my.integrate
        ### will use 1 sigma but cycle through thetas defined above
        sapply(X = theta, FUN = my.integrate, fun = integrand, sigma = sigma)
      } else {
        # use logistic approximation
        k = 16 * sqrt(3) / (15 * pi)
        SUMMER::expit(theta / sqrt(1 + k^2 * sigma^2))
      }
    }
  }
  
  # correcting post marginalisation if binomial
  ## do we want to use approximation or numerical integration for marginalisation
  if (inla_mod$family == 'binomial') {
    if(use.approximation){
      message("\n Binomial marginalisation found using approximation", appendLF = F)
      theta2 <- apply(X = cbind(tau, theta), 1, FUN = marginalise.logit.binomial, use.approximation = TRUE)
    } else{
      message("\n Binomial marginalisation found using numerical integration", appendLF = F)
      theta2 <- apply(X = cbind(tau, theta), 1, FUN = marginalise.logit.binomial, use.approximation = FALSE)
    }
    
  } else {
    theta2 <- SUMMER::expit(t(theta))
  }
  colnames(theta2) <- paste0('theta:', 1:ncol(theta2))
  
  # Suppress summarise info
  options(dplyr.summarise.inform = FALSE)
  
  # add additional columns
  ## proportions
  if(!is.null(strata.weights)){
    if(!is.null(inla_mod$Amat)){ # spatial
      if(sum(c('region', 'urban', 'period', 'prop') %in% colnames(strata.weights)) == 4){ # are correct columns supplied
        data3 <- 
          left_join(data2[,c('region', 'urban', 'strata', 'age', 'period', 'cohort')],
                    strata.weights, 
                    by = c('region', 'urban', 'period'))
        message("\n Subnational weights supplied", appendLF = F)
      } else {
        stop('Subnational weights supplied need to include: region, urban, period and prop column headings')
      }
    } else { # national
      if(sum(c('urban', 'period', 'prop') %in% colnames(strata.weights)) == 3){ # are correct columns supplied
        # not spatial
        data3 <- 
          left_join(data2[,c('region', 'urban', 'strata', 'age', 'period', 'cohort')], 
                    strata.weights, 
                    by = c('urban', 'period'))
        message("\n National weights included", appendLF = F)
      } else {
        stop('National weights supplied need to include: urban, period and prop column headings')
      }
    }
  } else {
    data3 <- 
      data2[,c('region', 'urban', 'strata', 'age', 'period', 'cohort')] %>% 
      mutate(prop = 0)
    message("\n No strata weights supplied. Set all weights to zero", appendLF = F)
  }
  ## labels from numbe of months in each age group
  labels <- as.character(unique(data2$age))
  ns <- rep(1, length(labels))
  for (i in 1:length(labels)) {
    if (labels[i] == "0") {
      ns[i] <- 1
      next
    }
    tmp <- as.numeric(strsplit(as.character(labels[i]), "-")[[1]])
    ns[i] <- tmp[2] - tmp[1] + 1
  }
  data3$ns <- ns[data3$age %>% as.numeric]
  
  # generic part of results
  ## true for overall and stratified
  ## true for national and subnational
  outTemp <- 
    cbind(data3, theta2) %>%
    # mutate(across(starts_with('theta:'), SUMMER::expit)) %>%
    group_by(region, urban, strata, age, period, prop, ns) %>%
    summarise(across(starts_with('theta:'), mean)) %>%
    ungroup %>%
    group_by(region, strata, period, prop) %>%
    # inside the product for survival
    # mutate(across(starts_with('theta:'), ~ ((1 - .x)^(ns)))) %>%
    mutate(across(starts_with('theta:'), ~ ((1 - .x)^ns))) %>%
    # take product over all age groups
    ungroup %>% group_by(region, strata, period, prop) %>%
    summarise(across(starts_with('theta:'), prod)) %>%
    # subract the survial probability from one to get mortality
    mutate(across(starts_with('theta:'), ~ 1 - .x)) %>%
    ungroup
  
  # stratified results
  ## by region and strata without proportions included
  out1 <- 
    outTemp %>%
    mutate(select(., starts_with('theta:')) %>% 
             apply(., 1, my.summary, lowerCI = lowerCI, upperCI = upperCI) %>% 
             lapply(., data.frame) %>% 
             do.call(rbind, .)) %>%
    select(-starts_with('theta:')) %>%
    as.data.frame()
  
  # overall results
  ## proportionally summed over strata
  out2 <-
    outTemp %>%
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
    select(-starts_with('theta:')) %>%
    as.data.frame()
  
  out <- list(stratified = out1, overall = out2)
  
  if(save.u5m.draws == TRUE){
    out$u5m.draws <- outTemp
  }
  
  end2 <- proc.time() - start2
  
  message("\nTime taken for posterior sampling:", 
          "\n Sampling (s):",  end1[3] %>% as.numeric %>% round(),
          "\n Cleaning (s):",  end2[3] %>% as.numeric %>% round(),
          "\n Total    (s):",  (end1[3] + end2[3]) %>% as.numeric %>% round(),
          "\n----------------------------------",
          appendLF = FALSE)
  
  return(out)
  
}

# function for summary
## quantiles, mean and varaince
my.summary <- function(x, lowerCI, upperCI) {
  
  # x = theta1[,1]; lowerCI = lowerCI; upperCI = upperCI
  
  qntl <- quantile(x, probs = c(lowerCI, 0.5, upperCI))
  data.frame(mean = mean(x), variance = var(x), lower = qntl[1], median = qntl[2], upper = qntl[3])
  
}

# map plots for APC results
my.map.plot <- function(data, variables, value = NULL, geo, by.data, by.geo, size = 0.5,  
                        border = "gray20", ncol = NULL, legend.label = NULL, 
                        color.option = c('A', 'B', 'C', 'D', 'E')[4], per1000 = FALSE,
                        direction = 1, removetab = FALSE, ylim = NULL) {
  
  # data =
  #   apcRes$stratified %>%
  #   filter(strata == 'rural') %>%
  #   mutate(region = region %>% as.character %>% tolower) %>%
  #   as.data.frame
  # variables = 'period';
  # by.data = 'region';
  # geo = geoPoly;
  # value = 'mean';
  # per1000 = removetab = FALSE;
  # clean = TRUE
  # by.geo = 'REGNAME'
  # ncol = 10;
  # direction = - 1;
  # legend.label = 'U5MR'
  # size = 0.5;
  # border = "gray20";
  # ylim = NULL;
  
  # select only variable, by data and value and reshape
  datTemp <- 
    data[, c(variables, by.data, value)] %>%
    reshape2::melt()
  
  # add in the correct label and scaling for variable and value
  datTemp$variable <- datTemp[,variables]
  if (per1000) {
    datTemp$value <- datTemp$value * 1000
  }
  
  # organise spatial locations
  has.coord <- !is.na(sp::proj4string(geo)) # coordinates for ggplot?
  geo2 <- ggplot2::fortify(geo, region = by.geo) # make a data frame from spatial object
  
  # combine geo data and model results
  data2 <- 
    merge(geo2, datTemp, by = 'id', by.y = by.data) %>% # need sort = FALSE here to stop polygon ordering being messed up
    arrange(order) # arrange by the fortifies objects order
  
  # make the plot
  ## need data and polygon
  p <- 
    ggplot2::ggplot(data2) +
    ggplot2::geom_polygon(aes(x = long, y = lat, group = group, fill = value), color = border, size = size)
  
  ## whether to show the tab label
  ## if only one variable then do not
  if (length(unique(data$variable)) > 1 || removetab == FALSE) {
    if (is.null(ncol)) {
      p <- p +  ggplot2::facet_wrap(~variable)
    } else {
      p <- p + ggplot2::facet_wrap(~variable, ncol = ncol)
    }
  }
  
  ## range of y values to be plotted
  if (!is.null(ylim)) {
    p <- p + ggplot2::scale_fill_viridis_c(legend.label, lim = ylim, direction = direction, option = color.option)
  } else {
    p <- p + ggplot2::scale_fill_viridis_c(legend.label, direction = direction, option = color.option)
  }
  
  ## add a coorinate map if avaliable
  if (has.coord){
    p <- 
      p + 
      ggplot2::coord_map()
  }
  
  return(p)
  
}

# theme for plots
my.theme<-function(...){
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black'),
        # legend.title=element_blank(),
        legend.text.align=0,
        legend.key=element_rect(fill=NA),
        ...)
}

# theme for map plots
my.map.theme <- function(...){
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 legend.text.align=0,
                 legend.key=element_rect(fill=NA),
                 panel.background=element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 ...)
}