library(data.table)
library(doBy)
library(dplyr)
library(sf)
library(spdep)

setwd('C:/fhi/excess_mortality/manuscript/nordics/PH/revision1/sw_submit/')
source('code/util.R')

req.nations <- c('DK', 'FI', 'SE', 
                 'BE', 'BG', 'ES', 'HU', 'PT', 'SK')

## model specification for INLA model
formula = 
  deaths ~ 1 + offset(log(population)) + 
  f(id.year, model='rw1') + 
  f(id.time, model='rw1', hyper=list(theta = list(prior="pc.prec", param=c(1, 0.01))), constr = TRUE, scale.model = TRUE, cyclic = TRUE) +
  f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, constr = TRUE, 
    hyper = list(theta1 = list('PCprior', c(1, 0.01)), theta2 = list('PCprior', c(0.5, 0.5)))) +
  f(id.space.year, model='iid')

control.family=inla.set.control.family.default()



### load all cause mortality dataset
d <- readRDS('data/Input_AllCauseMortality.RDS') |>
  filter(substring(geo.code, 1, 2) %in% req.nations) |> ## limit to 3 Nordic countries and 2 select European
  mutate(year=year(date), month=month(date)) 

d <- filter(d, !(substring(geo.code, 1, 2) =='DK' & date < '2007-01-01')) ## Remove years for whihc we have no data for Denmark.


  
# load shape file and limit to required countries. Needed by the INLA model
shp = read_sf("data/NUTS_RG_10M_2021_3035/")
shp = subset(shp, CNTR_CODE %in% req.nations) 
shp = subset(shp, LEVL_CODE == 2) 

shp = subset(shp, NUTS_ID %in% unique(d$geo.code)) ## limit to select countries
shp = shp[order(shp$NUTS_ID),]

# save shape file
W.nb <- poly2nb(shp)
nb2INLA("W.adj", W.nb) 


#### Build models ####
req.loc='SE' ## specify country we want to build models for


## Build INLA model. Takes <10 minutes. doINLA() filters out national data and only uses subregions
spatial.est <- splitBy(~age.group+sex, d[substring(geo.code, 1, 2) == req.loc]) |>
  lapply(doINLA) |>
  rbindlist() |>
  postProcessINLA(pop=d[substring(geo.code, 1, 2) == req.loc, 
                        j=.(pop=mean(population)), by=.(geo.code, geo.name, age.group, sex, year)]) ## reuse the dataset for population (annual=mean monthly(identical))


## Build GAMM model. Takes <10 minutes
gamm.est <- splitBy(~age.group+sex, d[geo.code == req.loc]) |>
  lapply(fitBayesianGAMM)

gamm.est <- list(
  lapply(gamm.est, function(x) x[[1]]) |> rbindlist(),
  lapply(gamm.est, function(x) x[[2]]) |> rbindlist()
)


## Build timeseries models. Takes <5 minutes
ts.est <- doAgeSex(d[geo.code == req.loc])



## stack estimates from the 3 models and calculate model average, monthly
ests.month <- rbind(mutate(gamm.est[[1]], q=as.numeric(q), model = 'GAMM') |> select(!c(pop, value)),
                    mutate(ts.est[[1]], month.date=as.Date(month.date), model = 'TS-combo')|> rename(date=month.date) |>
                      reshape2::melt(id.vars=c('model', 'geo.name', 'age.group', 'sex', 'date'), variable='q', variable.factor=F, value.name='rate'),
                    filter(spatial.est[[1]], geo.name %in% unique(gamm.est[[1]]$geo.name)) |> 
                      melt(id.vars=c('geo.name', 'age.group', 'sex', 'date'), variable='q', variable.factor=F, value.name='rate') |> 
                      mutate(q=as.numeric(q), model = 'INLA'))
ests.month <- rbind(ests.month, ests.month[, j=.(model='ENS_MEAN', rate=mean(rate, na.rm=T)), by=.(geo.name, age.group, sex, date, q)]) |>
  mutate(age.group=factor(age.group, levels=c('TOTAL', 'Y_LT60', 'Y60-69', 'Y70-79', 'Y_GE80')))

# stack estimates from the 3 models and calculate model average, annual
ests.yr <- rbind(mutate(gamm.est[[2]], q=as.numeric(q), model = 'GAMM') |> select(!c(pop, value)),
                 mutate(ts.est[[2]], model = 'TS-combo') |> 
                   reshape2::melt(id.vars=c('model', 'geo.name', 'age.group', 'sex', 'year'), variable='q', variable.factor=F, value.name='rate'),
                 filter(spatial.est[[2]], geo.name %in% unique(gamm.est[[1]]$geo.name)) |> 
                   melt(id.vars=c('geo.name', 'age.group', 'sex', 'year'), variable='q', variable.factor=F, value.name='rate') |> 
                   mutate(q=as.numeric(q), model = 'INLA')) 
ests.yr <- rbind(ests.yr, ests.yr[, j=.(model='ENS_MEAN', rate=mean(rate, na.rm=T)), by=.(geo.name, age.group, sex, year, q)]) |>
  mutate(age.group=factor(age.group, levels=c('TOTAL', 'Y_LT60', 'Y60-69', 'Y70-79', 'Y_GE80')))



#### Inspect estimates #####
require(ggplot2)
require(ggthemes)

theme_set(theme_clean() + 
            theme(legend.position = 'bottom', panel.background = element_rect(fill=NULL), 
                  panel.grid.major.x = element_line(linetype='dotted', color='darkgrey'),
                  axis.title = element_text(size=12, family='serif'), axis.text = element_text(size=10, family='serif'),
                  strip.text = element_text(size=14, family='serif')))


req.age='TOTAL' # all ages combined
req.sex='Total' # both sexes combined
req.colors <- c('black',  'darkblue', 'red', 'darkgreen')#  



## Monthly estimates per each model, by specified age and sex
filter(ests.month, age.group==req.age, sex == req.sex, q %in% c(0.025, .25, .5, .75, .975)) |>
  dcast(model+date~q, value.var = 'rate') |>
  ggplot(aes(x=date)) +
  geom_line(aes( y=`0.5`, color=model)) +
  geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`, fill=model), alpha=.4) +
  geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`, fill=model), alpha=.2) +
  geom_point(data=d[geo.code == req.loc & age.group==req.age & sex == req.sex & date > '2019-12-31'],
             aes(x=date, y=rate), shape=1, size=.9) +
  scale_color_manual(values=req.colors) +
  scale_fill_manual(values=req.colors) +
  facet_wrap(~model, ncol=1) +
  labs(x='Month', y='Mortality rate (per 100,000)', color='Model', fill='Model') 

## Monthly estimates per each model, by specified sex - age-stratified
filter(ests.month, sex == req.sex, q %in% c(0.025, .25, .5, .75, .975)) |>
  dcast(model+date+age.group~q, value.var = 'rate') |>
  ggplot(aes(x=date)) +
  geom_line(aes( y=`0.5`, color=model)) +
  geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`, fill=model), alpha=.4) +
  geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`, fill=model), alpha=.2) +
  geom_point(data=d[geo.code == req.loc & sex == req.sex & date > '2019-12-31'],
             aes(x=date, y=rate), shape=1, size=.9) +
  scale_color_manual(values=req.colors) +
  scale_fill_manual(values=req.colors) +
  facet_grid(age.group~model, scales='free_y') +
  labs(x='Month', y='Mortality rate (per 100,000)', color='Model', fill='Model')


## Annual estimates per each model, by specified sex and age
filter(ests.yr, age.group==req.age, sex == req.sex, q %in% c(0.025, .25, .5, .75, .975)) |>
  dcast(model+year~q, value.var = 'rate') |>
  ggplot() +
  geom_crossbar(aes(x=year, y=`0.5`, ymin=`0.025`, ymax=`0.975`, fill=model),
                width=.1, position=position_dodge(width=.5), alpha=.2, color=NA) +
  geom_crossbar(aes(x=year, y=`0.5`, ymin=`0.25`, ymax=`0.75`, fill=model, color=model),  
                width=.1, position=position_dodge(width=.5), alpha=.4, fatten=1.5) +
  geom_point(data=d[geo.code == req.loc & age.group==req.age & sex == req.sex & date > '2019-12-31',
                    j=.(rate=sum(deaths)*1E5/mean(population)), by=.(age.group, sex, year)],
             aes(x=year, y=rate), shape=1, size=3) +
  scale_x_continuous(breaks = 2020:2023) +
  scale_fill_manual(values=req.colors) +
  scale_color_manual(values=req.colors) +
  labs(x='Year', y='Mortality rate (per 100,000)', fill='Model') + guides(color='none') 


## Annual estimates per each model, by specified sex - age-stratified
filter(ests.yr, sex == req.sex, q %in% c(0.025, .25, .5, .75, .975)) |>
  dcast(model+age.group+year~q, value.var = 'rate') |>
  ggplot() +
  geom_crossbar(aes(x=year, y=`0.5`, ymin=`0.025`, ymax=`0.975`, fill=model),
                width=.1, position=position_dodge(width=.5), alpha=.2, color=NA) +
  geom_crossbar(aes(x=year, y=`0.5`, ymin=`0.25`, ymax=`0.75`, fill=model, color=model),  
                width=.1, position=position_dodge(width=.5), alpha=.4, fatten=1.5) +
  geom_point(data=d[geo.code == req.loc &  sex == req.sex & date > '2019-12-31',
                    j=.(rate=sum(deaths)*1E5/mean(population)), by=.(age.group, sex, year)],
             aes(x=year, y=rate), shape=1, size=1.5) +
  scale_x_continuous(breaks = 2020:2023) +
  scale_fill_manual(values=req.colors) +
  scale_color_manual(values=req.colors) +
  facet_wrap(~age.group, ncol=3, scales='free_y') +
  labs(x='Year', y='Mortality rate (per 100,000)', fill='Model') + guides(color='none') 

rm(list=ls())
gc()
