library(rLakeAnalyzer)
library(tidyverse)
library(lubridate)
library(reshape2)
library(scales)
library(gridExtra)
library(zoo)


# define functions

time_avg <- function(wt, period = "1 week"){
  # average water temperature data over a given time period
  # inputs:
  #  - wt: water temperature dataframe (must have datetime column called 'datetime')
  #  - period: time period to aggregate by (eg. '1 month' or '1 week')
  # output:
  #  - wt_avg: water temp data frame time averaged over the specified period
  time_bins = cut(wt$datetime, period)
  wt_avg <- aggregate(.~time_bins, mean, data=wt) %>%
    mutate(datetime=ymd(time_bins)) %>%
    select(!time_bins)
  return(wt_avg)
}


# function to create 'boxes' which are used for Kz calcs

create_layers <- function(bathy, interval=1.0, Zmax=NULL){
  # creates evenly spaced layers vertically in the water column
  # and outputs their key geometric properties
  # inputs: 
  # bathy: bathymetry of lake (data.frame of depth and area for given lake.
  # eg. from load.bathy, must be labelled 'depths' and 'areas'
  # interval: vertical height of layers, defaults to 1.0m
  # Zmax: max depth of lake. If not specified, will use max depth in bathys
  # output: data.frame?
  # with layer numbers, top and bottom layer depths,
  # areas of upper boundary of each layer and layer volumes
  
  # if depth not specified, use max depth in bathymetric data
  if(is.null(Zmax)){
    Zmax <- max(bathy$depths)
  }
  # upper and lower depths of each layer
  upperZ = seq(0.0, Zmax-interval, by = interval)
  lowerZ = seq(interval, Zmax, by = interval)
  # numbers for each layer (uppermost layer is = 1)
  numbers = 1:length(upperZ)
  
  # calc upper boundary area of each layer:
  # using linear interpolation between areas
  # assuming cone shape, and so add the zero area point at the bottom of the lake
  upperAreas = stats::approx(c(bathy$depths, Zmax), c(bathy$areas, 0.0), upperZ)$y
  lowerAreas = stats::approx(c(bathy$depths, Zmax), c(bathy$areas, 0.0), lowerZ)$y
  # calc volume of each layer
  # currently just average upper and lower area and then multiply by layer height
  # BUT actually looks like bathymetry data sometimes comes with layer volumes! 
  # so could add this as an optional input. Might be interesting as well to see
  # how well my volume calc method approximates the volumes given for blelham
  volumes = interval * (upperAreas + lowerAreas) / 2.0
  
  return(data.frame(layer_num=numbers, upperZ=upperZ, lowerZ=lowerZ, 
                    upperAreas=upperAreas, lowerAreas=lowerAreas, volume=volumes))
}

layer_temps <- function(wtr_temps, layers, bathy){
  # get avg temperatures of evenly spaced layers in the timeseries
  # inputs:
  # wtr_temps - data.frame in standard Rlakeanalyser timeseries format 
  # layers - data.frame from create_layers() function
  #   must include upperZ, lowerZ columns
  # bathy - from load.bathy. Could technically use data from layers though
  # outputs: data.frame with avg layer temps through entire timeseries
  # 
  layer_tmps = tibble(datetime = wtr_temps$datetime)
  for (i in 1:nrow(layers)){
    layer_tmp = ts.layer.temperature(wtr_temps, layers$upperZ[i],
                                     layers$lowerZ[i], bathy)
    names(layer_tmp)[2] <- paste0('wtr_',
                                  as.character(layers$upperZ[i]))
    layer_tmps <- left_join(layer_tmps, layer_tmp, by='datetime')
    
  }
  return(layer_tmps)
  
  
}


tempGradient <- function(wtr_temps){
  # approximate temperature gradients with respect to time
  # this function outputs a single gradient for each time step,
  # calculated by looking at the difference between it and the previous
  # inputs:
  # wtr_temps - data.frame in standard Rlakeanalyser timseries format 
  # must have a 'datetime' column and 'wtr_x' columns for the temperature at depth x
  # outputs: tibble with columns for datetime, temp_gradient_x
  wtr_temps$datetime = as.Date(wtr_temps$datetime)
  
  wtr_time_grads <- wtr_temps %>%
    # find time step between measurements
    mutate(date_diff = difftime(datetime, lag(datetime), units='secs')) %>%
    # for all wtr temps (numerics) find diff and divide by date_diff
    mutate(across(
      where(is.numeric), (function(x) x - lag(x))) / as.numeric(date_diff)) %>%
    # remove redundant date_diff variable
    select(!date_diff) %>%
    # remove NAs (which arise as first variable has no previous value for grad calculation)
    drop_na() %>%
    # put into long format as I think this will make things easier!
    pivot_longer(cols = where(is.numeric),
                 names_to = "upperZ",
                 values_to = "gradient") %>%
    separate(upperZ, c(NA, 'upperZ'), sep = '_') %>%
    # extract numerical value from depth
    mutate(upperZ=as.numeric(upperZ))
  return(wtr_time_grads)
}

tempGradientSingle <- function(wtr_temps, dateStart=NULL, dateStop=NULL){
  # approximate temperature gradients with respect to time over a date range
  # this function outputs a SINGLE time-averaged gradient for each depth! 
  # inputs:
  # wtr_temps - data.frame in standard Rlakeanalyser timseries format 
  # must have a 'datetime' column and 'wtr_x' columns for the temperature at depth x
  # dateStart, dateStop - optional start and end dates for analysis
  #   if not specified, then function will use entire daterange
  # outputs: tibble with columns for depth, gradient, and adjusted R^2 values
  # 
  # need datetime in POSIXct class in order for correct gradient calculation
  wtr_temps$datetime = as.Date(wtr_temps$datetime)
  
  if(is.null(dateStart)){
    dateStart = min(wtr_temps$datetime)
  }
  if(is.null(dateStop)){
    dateStop = max(wtr_temps$datetime)
  }
  # select data in specified range
  wtr.subset = subset(wtr_temps, datetime >= dateStart & datetime <= dateStop)
  
  # calc gradient of temp/time data using linear regression
  wtr.grad <- wtr.subset %>%
    summarise(
      across(where(is.numeric), ~lm(.x ~ datetime)$coefficients[2], 
             .names="grad-{.col}")) %>%
    pivot_longer(cols = starts_with('grad'), names_to = 'depth', 
                 values_to = 'gradient') %>%
    separate(col = depth, sep = '_', into = c(NA, 'depth')) %>%
    mutate(depth = as.numeric(depth))
  
  # calc adj R^2 values for linear regression models:
  wtr.Rsq <- wtr.subset %>%
    summarise(
      across(where(is.numeric), ~summary(lm(.x ~ datetime))$adj.r.squared, 
             .names="adjR^2-{.col}")) %>%
    pivot_longer(cols = starts_with('adj'), names_to = 'depth', 
                 values_to = 'adj_R^2') %>%
    separate(col = depth, sep = '_', int = c(NA, 'depth')) %>%
    mutate(depth = as.numeric(depth))
  
  return(left_join(wtr.grad, wtr.Rsq, by='depth'))
}

vertTempGrad <- function(layer_tmps, layers){
  # calculate vertical temperature gradients between layers
  # inputs:
  # layer_tmps - avg temperatures for each layer, from layer_temps() function
  # layers - layer data created from create_layers()
  # outputs: data frame with columns for datetime, G (vert temp grad), 
  #                                      and upperZ (depth of upper boundary)
  #          
  
  # initialise as tibble
  G = tibble(datetime = as.Date(layer_tmps$datetime))
  
  for (i in 2:nrow(layers)){
    Z_above = layers$upperZ[i-1]
    Z = layers$upperZ[i]
    deltaZ = Z - Z_above
    temps_above = select(layer_tmps,
                         ends_with(paste0('_', as.character(Z_above))))
    temps_Z = select(layer_tmps,
                     ends_with(paste0('_', as.character(Z))))
    G_name = paste0('G_', as.character(Z))
    
    G[G_name] = (temps_above - temps_Z)/deltaZ
  }
  G <- G %>% pivot_longer(where(is.numeric), 
                          names_to = c(NA, 'upperZ'), names_sep = '_', 
                          values_to = 'G') %>%
    mutate(upperZ=as.numeric(upperZ))
  return(G)
  
}

calcKzOLD <- function(wtr, bathy, deltaZ=1.0, Zmax=NULL){
  # calc Kzs of water temperature timeseries
  # inputs:
  # wtr - data.frame of water temps in standard Rlakeanalyser timeseries format
  # bathy - bathymetric data from load.bathy
  # deltaZ - depth of layers to evaluate Kz over. Defaults to 1.0 m
  # maxZ - maximum depth of lake. 
  #        If not specified, will use max depth from bathymetry
  # output: timeseries of Kz values for each layer.
  
  # create layers:
  layers = create_layers(bathy, interval = deltaZ, Zmax=Zmax)
  # get layer temps
  layer.temps = layer_temps(wtr, layers, bathy)
  # find temperature time gradients:
  time.grads = tempGradient(layer.temps)
  # find vertical gradients, G
  G = vertTempGrad(layer.temps, layers)
  # intialise Kz data frame
  Kz = tibble(datetime = layer.temps$datetime)
  m = max(layers$layer_num)
  for (i in 2:m){
    Z = layers$upperZ[i]
    Kzname = paste0('Kz_', as.character(Z))
    dTdt = subset(time.grads, time.grads$depth >= Z)$gradient
    vols = subset(layers, layers$upperZ >= Z)$volume
    Gz = select(G, ends_with(paste0('_', as.character(Z))))
    Kz[Kzname] = sum(dTdt * vols) / (Gz * layers$upperAreas[i])
  }
  return(Kz)
}

calcKz <- function(wtr, bathy, deltaZ=1.0, Zmax=NULL){
  # calc Kzs of water temperature timeseries
  # inputs:
  # wtr - data.frame of water temps in standard Rlakeanalyser timeseries format
  # bathy - bathymetric data from load.bathy
  # deltaZ - depth of layers to evaluate Kz over. Defaults to 1.0 m
  # maxZ - maximum depth of lake. 
  #        If not specified, will use max depth from bathymetry
  # output: timeseries of Kz values for each layer.
  
  # create layers:
  layers = create_layers(bathy, interval = deltaZ, Zmax=Zmax)
  # get layer temps
  layer.temps = layer_temps(wtr, layers, bathy)
  # find temperature time gradients:
  time.grads = tempGradient(layer.temps)
  # find vertical gradients, G
  G = vertTempGrad(layer.temps, layers)
  
  Kz_DF <- layers %>% left_join(time.grads) %>%
    left_join(G) %>% drop_na() %>%
    # sort by descending upperZ
    arrange(desc(upperZ)) %>%
    group_by(datetime) %>%
    # calc sums of vertical gradients with vols 
    # (top part of fraction in equation)
    mutate(sumgradXvol = cumsum(gradient * volume)) %>%
    group_by(datetime, upperZ) %>%
    # divide by lower part of fraction in equation
    summarise(Kz = sumgradXvol/(upperAreas*G))
  
  return(Kz_DF)
}

# note, the following function is useful for calcKzOLD
# but isn't needed for current calcKz method
pivot_longer_wtrDF <- function(wtrDF, colStart='Kz_1', 
                               colEnd='Kz_13', sep='_'){
  # pivots the input 'wide' tibble wtrDF into 'long' format
  # inputs:
  # - wtrDF: water data frame, in wide format with variables measured
  #          at multiple depths (eg. Kz_1, Kz_2 ...)
  #          depth (m) must be included in variable names separated by sep (eg. '_')
  # - colStart: name of first column to widen over
  # - colEnd: name of last column to pivot over
  # - sep: text separator used to separate depth from var name 
  # note that all columns between colStart and colEnd will be pivoted
  
  wt.tidy <- wtrDF %>%
    # create separate rows for each depth
    # forced colStart and colEnd to be characters because it was thinking 
    pivot_longer(cols = as.character(colStart):as.character(colEnd),
                 names_to = "depth",
                 values_to = str_split_fixed(colStart, pattern=sep, 2)[1,1]) %>%
    separate(depth, c(NA, 'depth'), sep = sep) %>%
    # extract numerical value from depth
    mutate(depth=as.numeric(depth))
  return(wt.tidy)
}

upsample_ts <- function(wt, period = '1 day', start=NULL, keepOrig=FALSE){
  # upsample lake timeseries data to a given period
  # this will actually downsample too! 
  #  - wt: input dataframe (must have datetime column called 'datetime')
  #  - period: time period to upsample to (eg. '1 month' or '1 week')
  #  - start: first date in timeseries, to allow for sync with other timeseries 
  #           but, if left as NULL then will be set to first datetime in wt
  #  - keepOrig: whether to keep original values that don't fit into the new 
  #              timeseries dates
  # output: 
  #  - wt_up: upsampled timeseries data
  
  # define start and end of dates
  if(is.null(start)){
    start = wt$datetime[1]
  }
  # force start to be datetime:
  start = as.Date(start)
  # check start date is valid
  if(start < wt$datetime[1]){
    print('Given start date is before first observation!')
  }
  end = as.Date(tail(wt$datetime, 1))
  # create sequence of dates with upsampling period
  date.indexes = seq(
    from = start,
    to = end,
    by = period
  )
  date.indexes = tibble(datetime = date.indexes)
  head(date.indexes)
  # force datetime of wt to be posixct
  wt$datetime = as.Date(wt$datetime)
  # join upsampled date sequence with data
  # this produces NAs for all new dates
  wt_up <- full_join(date.indexes, wt, by='datetime') %>%
    # use na.approx to fill NAs with interpolated values
    # for all numeric columns 
    # NOTE: this works, but interpolation may not be the best approach
    # could use spline fitting (ie. non-linear)
    # (eg. just fill missing values with previous data value)
    mutate(across(
      where(is.numeric), (function(x) na.approx(x, datetime)) )) %>%
    arrange(datetime)
  if(!keepOrig){
    # use left join to remove entries with dates not in date.indexes
    wt_up <- left_join(date.indexes, wt_up, by='datetime')
  }
  return(wt_up)
}

ts.water.density <- function(wtr, sal = wtr * 0){
  # function to calculate water density of a timeseries of water temps
  # inputs: 
  #   wtr - water temperature data frame in standard lake analyser format
  #   sal - salinity, defaults to 0 (currently not implemented!)
  # outputs: timeseries of densities at different depths in long format
  wtr.long <- pivot_longer_wtrDF(wtr, colStart = colnames(wtr)[2], 
                                 colEnd = tail(colnames(wtr), n=1), sep='_') %>%
    mutate(density=water.density(wtr))#, sal = sal))
  
  return(wtr.long)
}

ts.heat.flux <- function(KzDF, wtr, bathy, 
                         deltaZ=1.0, Zmax=NULL, cp = 4184, sal = wtr*0){
  # calculate heat fluxes at top of each layer from time series of Kz values, 
  # water temp readings,and bathymetry
  # inputs:
  #   Kz - Kz data from calcKz()
  #   wtr - water temp dataframe in standard wide format
  #   bathy - bathymetric data from load.bathy
  #   deltaZ - depth of layers to evaluate Kz over. Defaults to 1.0 m
  #   Zmax - maximum depth of lake. 
  #        If not specified, will use max depth from bathymetry
  #   cp - specific heat capacity of water at const. pressure 
  #          (default=4184 J/kg/K)
  #   sal - salinity data, defaults to 0 (currently not implemented!)
  # outputs: heat flux data frame time series
  
  
  # create layers:
  layers = create_layers(bathy, interval = deltaZ, Zmax=Zmax)
  # get layer temps
  layer.temps = layer_temps(wtr, layers, bathy)
  # find vertical gradients, G
  G = vertTempGrad(layer.temps, layers)
  # calc densities
  density = ts.water.density(layer.temps) %>%#, sal) %>%
    rename(upperZ=depth)
  
  Kz.d.G.j <- left_join(KzDF, density) %>%
    left_join(G) %>%
    mutate(j = Kz * density * G * cp)
  
  return(Kz.d.G.j)
}