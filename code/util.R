library(data.table)
library(doBy)
library(dplyr)

library(rstanarm) # GAMM implementation

library(INLA) # for INLA model
library(sf)
library(spdep)

library(fable) # for TS-Combo model
library(fabletools)
library(feasts)


## Fits GAMM model
fitBayesianGAMM <- function(d, train.end='2020-01-01'){
  train <- subset(d, date < train.end)
  test <- subset(d, date >= train.end)
  
  train <- train[date >= '2000-01-01']
  
  
  model <- rstanarm::stan_gamm4(
    deaths ~ offset(log(population)) + # age-specific population offset
      s(year, k = length(unique(train$year))) + # annual trend
      s(month, bs = "cc", k = length(unique(train$month))) , 
    data = train,
    ## Set family. Possibly neg_binomial_2(), poisson()or gaussian() 
    chains = 2, # increased from default =1. problems converging with fewer chains
    seed = 3542
  )
  
  
  ### Inspect model
  # summary(model)
  # plot_nonlinear(model)
  # pp_check(model, nreps=1e2, plotfun = "ppc_hist", seed=3542)
  # shinystan::launch_shinystan(model) ### fuller exploration
  
  
  ### Get posterior over projection period
  pred <- rstanarm::posterior_predict(
    model,
    newdata=test,
    offset=log(test$population)
  )
  
  pred.o <- melt(data.table(subset(test, select=c(age.group, sex, geo.name, year, date, population)), t(pred)), 
                 id.vars = c('age.group', 'sex', 'geo.name', 'year', 'date', 'population'))
  
  ### from posterior calculate estimates at the required quantile levels
  monthly <- pred.o[, j=.(value=sum(value), pop=sum(population, na.rm=T)), 
                    by=.(sex, age.group, date, variable)][, 
                                                          j=.(`0.5` = quantile(value, .5), 
                                                              `0.25` = quantile(value, .25), `0.75` = quantile(value, .75), 
                                                              `0.2` = quantile(value, .2), `0.8` = quantile(value, .8), 
                                                              `0.1` = quantile(value, .1), `0.9` = quantile(value, .9), 
                                                              `0.05` = quantile(value, .05), `0.95` = quantile(value, .95),
                                                              `0.025` = quantile(value, .025), `0.975` = quantile(value, .975)), 
                                                          by=.(date, sex, age.group, pop)] |>
    melt(id.vars=c('date', 'sex', 'age.group','pop'), variable.factor = F, variable.name = 'q') |>
    mutate(rate=value*1E5/pop,
           geo.name=unique(d$geo.name))
  
  
  
  #####
  ### aggregate posterior monthly to annual and calculate estimates at the required quantile levels
  yearly <- pred.o[, j=.(value=sum(value), pop=sum(population, na.rm=T)), 
                   by=.(sex, age.group, date, year, variable)][, 
                                                               j=.(value=sum(value), pop=mean(pop)),
                                                               by=.(sex, age.group, year, variable)][, 
                                                                                                     j=.(`0.5` = quantile(value, .5), 
                                                                                                         `0.25` = quantile(value, .25), `0.75` = quantile(value, .75), 
                                                                                                         `0.2` = quantile(value, .2), `0.8` = quantile(value, .8), 
                                                                                                         `0.1` = quantile(value, .1), `0.9` = quantile(value, .9), 
                                                                                                         `0.05` = quantile(value, .05), `0.95` = quantile(value, .95),
                                                                                                         `0.025` = quantile(value, .025), `0.975` = quantile(value, .975)), 
                                                                                                     by=.(year, sex, age.group, pop)]%>%
    melt(id.vars=c('year', 'sex', 'age.group','pop'), variable.factor = F, variable.name = 'q') |>
    mutate(rate=value*1E5/pop,
           geo.name=unique(d$geo.name))
  
  
  list(monthly, yearly) ## list of two objects returned
}



# Fits TS_Combo model
doAgeSex <- function(x, qs=c(.025, .05, .1, .2, .25, .5, .75, .8, .9, .95, .975), req.h = '4 years'){
  pop <- x[, j=.(pop=mean(population)), by=.(geo.name, age.group, sex, year) ]
  
  # limit to pre-pandemic years. 
  x <- orderBy(~date, x[date < '2020-01-01' & date >='2000-01-01'])
  x$month.date <- tsibble::yearmonth(gsub('-','',x$date))
  
  # fit model. note specification of key. This is important.
  # We treat each location+age+sex as separate time series
  # also getting a simple average of the 2 models
  ## See Chapter 11 of Hyndman et al FPP3: https://otexts.com/fpp3/hierarchical.html, for potential refinements. Forecast reconciliation etc. 
  fit <- x  |>
    as_tsibble(index=month.date, key= c(geo.name, age.group, sex)) |>
    model(arima = ARIMA(rate~pdq()+PDQ()+season(period=12)), 
          stl_ets = decomp_spec) |>
    mutate(combination = (arima + stl_ets)/2)
  
  
  
  ### Aggregation starts here. We first retain monthly estimates, but later also aggregate to year ####
  ### See https://otexts.com/fpp3/aggregates.html, which recommends using simulations. 
  fcast <- fit |> 
    select(!c(arima, stl_ets)) |> # retain only combination
    generate(times = 1000, h = req.h)  |>
    rename(rate=.sim) |>
    as_tibble() 
  
  fcast.agg <-  fcast |>
    group_by(geo.name, age.group, sex, month.date) |>
    summarise(rate = distributional::dist_sample(list(rate)))  |>  # Store as a distribution
    mutate(`0.5` = quantile(rate, .5), 
           `0.25` = quantile(rate, .25), `0.75` = quantile(rate, .75), 
           `0.2` = quantile(rate, .2), `0.8` = quantile(rate, .8), 
           `0.1` = quantile(rate, .1), `0.9` = quantile(rate, .9), 
           `0.05` = quantile(rate, .05), `0.95` = quantile(rate, .95),
           `0.025` = quantile(rate, .025), `0.975` = quantile(rate, .975)
    )  |>
    select(!rate)
  
  
  
  
  ### As above, but aggregating to YEAR ####
  ### here we are only collapsing across weeks in a year, so the denominator is mean pop
  fcast.yr.agg <- merge(mutate(fcast, year = year(month.date)), 
                        pop, 
                        by = c('geo.name', 'age.group', 'sex', 'year')) |>
    group_by(geo.name, age.group, sex, year, .rep) |>
    summarise(rate=sum(rate*pop/1E5)*1E5/mean(pop))|>
    group_by(geo.name, age.group, sex, year) |>
    summarise(rate = distributional::dist_sample(list(rate)))  |>  # Store as a distribution
    mutate(`0.5` = quantile(rate, .5), 
           `0.25` = quantile(rate, .25), `0.75` = quantile(rate, .75), 
           `0.2` = quantile(rate, .2), `0.8` = quantile(rate, .8), 
           `0.1` = quantile(rate, .1), `0.9` = quantile(rate, .9), 
           `0.05` = quantile(rate, .05), `0.95` = quantile(rate, .95),
           `0.025` = quantile(rate, .025), `0.975` = quantile(rate, .975)
    )  |>
    select(!rate)
  
  list(fcast.agg, fcast.yr.agg)
}



## TS_Combo helper: STL-ETS model definition
decomp_spec <- decomposition_model(
  STL(rate~trend()+season(period=12)),
  ETS(season_adjust ~ season("N"))
)

## TS_Combo helper: Get median, IQR and 95% intervals. Expand quantile vector to get additional estimates
getQuantiles <- function(f, req.q){
  lapply(req.q, function(y){
    mutate(f, q=y, rate=quantile(rate, y), .keep='unused')
  }) |>
    rbindlist()
}



### Build INLA model
doINLA <- function(d.sub){
  d.sub <- d.sub[geo.code %in% unique(shp$NUTS_ID)] ## this removes our national data and retain subnational only
  
  d.sub[is.na(deaths)]$deaths <- 0 # if no recorded deaths assume 0
  
  # Create indices
  d.sub <- mutate(d.sub, 
                  id.space = as.numeric(as.factor(geo.code)),
                  id.year = year - min(year) + 1,
                  id.time = month,
                  id.space.year = paste0(id.space, '_', id.year)) |>
    arrange(id.space)
  
  
  # subset data for 2020 and later
  truth <- d.sub[year>=2020]
  
  # set 2020 and later data to NA in training set. This specifies the estimates we want
  d.sub[year >= 2020]$deaths <- NA
  
  ret <- NULL
  
  # run model with try
  mod = try(inla(formula,
                 data=d.sub,
                 family="poisson", #poisson, nbinomial or see names(inla.models()$likelihood)  
                 verbose = F, 
                 control.family=control.family,
                 control.compute=list(config = TRUE, dic = TRUE, waic = TRUE), 
                 control.mode=list(restart=T),
                 num.threads = round(parallel::detectCores()*.8), 
                 control.predictor = list(link = 1)))
  
  ### Inspect model
  # summary(mod) ## overall summary
  # mod$summary.hyperpar ## hyperparameters
  # mod$summary.fixed
  # mod$summary.random
  
  if(!class(mod) == 'try-error'){ # if successful
    set.seed(3542)
    post.samples <- inla.posterior.sample(n = 1000, result = mod)
    predlist <- do.call(cbind,
                        lapply(post.samples,
                               function(x) exp(x$latent[startsWith(rownames(x$latent), "Pred")])))
    
    pois.samples <- apply(predlist, 2, function(y) rpois(n = length(y), lambda = y)) |>
      as.data.frame()
    
    # use if model was built with family = nbinomial. size = 1/dispersion is a hyperparam
    # pois.samples <- apply(predlist, 2, function(y) rnbinom(n = length(y), mu = y, size=mod$summary.hyperpar[1, '0.5quant'])) |>
    #   as.data.frame()
    
    ## add keys
    pois.samples$geo.code <- d.sub$geo.code
    pois.samples$year <- d.sub$year
    pois.samples$date <- d.sub$date
    
    ## only retain predictions (drop fit)
    pois.samples <- subset(pois.samples, year>=2020)
    
    
    pois.samples <- merge(pois.samples, 
                          subset(truth, select = -c(id.space, id.year, id.time, id.space.year)), 
                          by = c('year', 'date', 'geo.code'))  |>
      data.table()
    
    
    
    
    ret <- pois.samples
    rm(pois.samples, post.samples, predlist, mod)
    gc()
    gc()
    
  }
  
  return(ret) ## return full posterior
  
}

### Helper functions for INLA ####
postProcessINLA <- function(ests, pop, sample.size=1E3){
  inla.yr <- rbind(
    doCauseSub(ests, by.var = c('age.group', 'sex')) |> mutate(geo.name=unique(pop[nchar(geo.code) == 2]$geo.name)), ## National estimates by age and sex
    doCauseSub(ests, by.var = c('geo.name', 'age.group', 'sex')), ## subregion estimates by age and sex
    fill=T) |>
    select(-group.by)
  
  inla.yr <- merge(melt(inla.yr, 
                        id.vars=c('geo.name', 'age.group', 'sex', 'year'),
                        variable.name = 'q', value.name = 'expected'),
                   pop, 
                   by=c('geo.name', 'age.group', 'sex', 'year')) |>
    mutate(exp.rate = round(expected * 1E5/pop, 4))
  inla.yr <- dcast(inla.yr, geo.name+age.group+sex+year~q, value.var = 'exp.rate')
  
  
  
  
  ## monthly subnational
  inla.smry <- ests |> 
    rowwise(date) |> 
    summaryDists() |> 
    select(!paste0('V', 1:sample.size)) |> 
    data.table()
  
  # monthly national
  inla.month <- doCauseSub.month(ests, 1e3, c('age.group', 'sex')) |> 
    mutate(geo.name=unique(pop[nchar(geo.code) == 2]$geo.name))  |>
    select(-group.by)
  
  # stack subnational and national
  inla.month <- rbind(inla.month, select(inla.smry, names(inla.month))) |>
    melt(id.vars=c('geo.name', 'age.group', 'sex', 'date'),
         variable.name = 'q', value.name = 'expected') |>
    mutate(year=year(date)) |>
    merge(pop, by=c('geo.name', 'age.group', 'sex', 'year')) |>
    mutate(exp.rate = round(expected * 1E5/pop, 4))
  
  inla.month <- dcast(inla.month, geo.name+age.group+sex+year+date~q, value.var = 'exp.rate')
  inla.month$year = NULL
  
  return(list(inla.month, inla.yr))
}




doCauseSub <- function(post, sample.size=1e3, by.var=NULL){
  post <- splitBy(~year, post) |>
    lapply(function(x){
      splitBy(by.var, x) |> 
        lapply(groupCriterion, req.sample.size=sample.size) |>
        rbindlist(idcol='group')
    }) |>
    rbindlist(idcol='year') |>
    mutate(year=as.integer(year))
  
  
  post <- collapseVars(post, by.var) |>
    select(!group) |>
    mutate(group.by = paste(by.var, collapse='|'))
  
  
  by.var.2 <- c('year', 'group.by', by.var)
  post |>
    rowwise(all_of(by.var.2)) |> 
    summaryDists() |> 
    select(!paste0('V', 1:sample.size)) |>
    data.table()
}


doCauseSub.month <- function(post, sample.size=1e3, by.var=NULL){
  post <- splitBy(~date, post) |> 
    lapply(function(x){
      splitBy(by.var, x) |>
        lapply(groupCriterion, req.sample.size=sample.size) |>
        rbindlist(idcol='group')
    }) |>
    rbindlist(idcol='date')|>
    mutate(date = as.Date(date))
  
  post <- collapseVars(post, by.var) |>
    select(!group) |>
    mutate(group.by = paste(by.var, collapse='|'))
  
  
  by.var.2 <- c('date', 'group.by', by.var)
  post |> 
    rowwise(all_of(by.var.2)) |> 
    summaryDists() |> 
    select(!paste0('V', 1:sample.size)) |>
    data.table()
  
}


## summarize posteriors, get quantile estimates at required levels
summaryDists <- function(x){
  x |>
    mutate(`0.5` = median(c_across(V1:V1000)), 
           `0.25` = quantile(c_across(V1:V1000), probs= 0.25), 
           `0.75` = quantile(c_across(V1:V1000), probs= 0.75),
           `0.2` = quantile(c_across(V1:V1000), probs= 0.2), 
           `0.8` = quantile(c_across(V1:V1000), probs= 0.8),
           `0.1` = quantile(c_across(V1:V1000), probs= 0.1), 
           `0.9` = quantile(c_across(V1:V1000), probs= 0.9),
           `0.05` = quantile(c_across(V1:V1000), probs= 0.05), 
           `0.95` = quantile(c_across(V1:V1000), probs= 0.95),
           `0.025` = quantile(c_across(V1:V1000), probs= 0.025), 
           `0.975` = quantile(c_across(V1:V1000), probs= 0.975)
    )
}

collapseVars <- function(x, vars){
  group.names <- strsplit(unique(x$group), '\\|') |> unlist()
  group.names <- matrix(group.names, ncol=length(vars), byrow = T) |> data.table()
  
  setnames(group.names, vars)
  group.names$group=unique(x$group)
  
  merge(x, group.names, by='group')
}


groupCriterion <- function(x, req.sample.size=1e3){
  select(x, starts_with("V") & !matches('value')) |>
    apply(2, sum) |>
    matrix(ncol=req.sample.size) |>
    as.data.frame()
}

