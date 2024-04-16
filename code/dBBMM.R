# libraries and establishing API key -------------------------------------------------------

library(movegroup)
library(ggmap)

register_google(key = "###", write = TRUE) # rerun this periodically with new API key

# dusky dBBMM -------------------------------------------------------------

duskySaveDir <- paste0(owd, "/duskyAllRasters") #new folders for each species and 

### could potentially make this a for loop by making a list of each deployment (or year) and each species
### then it would filter through all combos, creating separate rasters for each 
dusky_dBBMM <- orstedPosMerged2 %>% 
  filter(common_name_e == "Dusky") %>% 
  mutate(est = as.POSIXct(est), 
         time = as.POSIXct(time),
         name = as.factor(name)) %>% 
  arrange(est) %>% 
  arrange(name)
str(dusky_dBBMM)

### When used together, the order of functions would be: movegroup, scaleraster, alignraster if required, plotraster.
### movegroup is nice because it removes individuals that do not have enough positions/detections
movegroup(
  data = dusky_dBBMM,
  ID = "name",
  Datetime = "est",
  Lat = "latitude",
  Lon = "longitude",
  savedir = duskySaveDir) #works with sandbar, sets temporary working directory for individual lever UDs

# Having run the movegroup function example:
scaleraster(path = duskySaveDir) 

# Weighted by number of positions per ID, fewer locations = lower Weighting value = higher final 
  # UD values after dividing by Weighting. This scales all IDs up to match the group max.
  # I think this is only necessary if there is imbalance in receiver locations but since I'm working with positions not really necessary

#scaleraster(path = tempdir(), weighting = Weighting)

#alginraster combines region-specific group-level UD into a single raster

stationLatsD2 <- stations %>% 
  filter(dep_period == "D1 target") %>% 
  select(latitude) %>% 
  dplyr::rename(lat = "latitude")

stationLonD2 <- stations %>% 
  filter(dep_period == "D1 target") %>% 
  select(longitude) %>% 
  dplyr::rename(lon = "longitude")

plotraster(
  x = paste0(duskySaveDir, "/Scaled/All_Rasters_Scaled_Weighted_UDScaled.asc"),
  mapsource = "google",
  maptype = "satellite",
  savedir = paste0(duskySaveDir, "/Plot"),
  xlatlon = paste0(duskySaveDir, "/Scaled/All_Rasters_Scaled_Weighted_LatLon.asc"),
  locationpoints = dusky_dBBMM |> dplyr::rename(lat = "latitude", lon = "longitude"),
  #receiverlats = subset(stations, stations$dep_period == "D2 target") |> 
    #dplyr::rename(lat = "latitude"),
  receiverlats = stationLatsD2$lat,
  receiverlons = stationLonD2$lon, 
  #receiverlons = subset(stations, stations$dep_period == "D2 target") |> 
    #dplyr::rename(lon = "longitude"),
  recpointsfill  = "white",
  recpointssize = 3,
  plotsubtitle = expression(paste("Space use by dusky sharks ", italic("Carcharhinus obscurus"))),
  pointsincontourssave = paste0(duskySaveDir, "/Scaled/pointsincontours.csv"))  

# sand tiger dBBMM --------------------------------------------------------

sandtigerSaveDir <- paste0(owd, "/sandtigerAllRasters") #new folders for each species and 

### could potentially make this a for loop by making a list of each deployment (or year) and each species
### then it would filter through all combos, creating separate rasters for each 
sandtiger_dBBMM <- orstedPosMerged2 %>% 
  filter(common_name_e == "Sand Tiger") %>% 
  mutate(est = as.POSIXct(est), 
         time = as.POSIXct(time),
         name = as.factor(name)) %>% 
  arrange(est) %>% 
  arrange(name)
str(sandtiger_dBBMM)

### When used together, the order of functions would be: movegroup, scaleraster, alignraster if required, plotraster.
### movegroup is nice because it removes individuals that do not have enough positions/detections
movegroup(
  data = sandtiger_dBBMM,
  ID = "name",
  Datetime = "est",
  Lat = "latitude",
  Lon = "longitude",
  savedir = sandtigerSaveDir) #works with sandbar, sets temporary working directory for individual lever UDs

# Having run the movegroup function example:
scaleraster(path = sandtigerSaveDir) #owd = output working directory established early on

# Weighted by number of positions per ID, fewer locations = lower Weighting value = higher final 
# UD values after dividing by Weighting. This scales all IDs up to match the group max.
# I think this is only necessary if there is imbalance in receiver locations but since I'm working with positions not really necessary

#scaleraster(path = tempdir(), weighting = Weighting)

#alginraster combines region-specific group-level UD into a single raster

stationLatsD2 <- stations %>% 
  filter(dep_period == "D1 target") %>% 
  select(latitude) %>% 
  dplyr::rename(lat = "latitude")

stationLonD2 <- stations %>% 
  filter(dep_period == "D1 target") %>% 
  select(longitude) %>% 
  dplyr::rename(lon = "longitude")

plotraster(
  x = paste0(sandtigerSaveDir, "/Scaled/All_Rasters_Scaled_Weighted_UDScaled.asc"),
  mapsource = "google",
  maptype = "satellite",
  savedir = paste0(sandtigerSaveDir, "/Plot"),
  xlatlon = paste0(sandtigerSaveDir, "/Scaled/All_Rasters_Scaled_Weighted_LatLon.asc"),
  locationpoints = sandtiger_dBBMM |> dplyr::rename(lat = "latitude", lon = "longitude"),
  #receiverlats = subset(stations, stations$dep_period == "D2 target") |> 
  #dplyr::rename(lat = "latitude"),
  receiverlats = stationLatsD2$lat,
  receiverlons = stationLonD2$lon, 
  #receiverlons = subset(stations, stations$dep_period == "D2 target") |> 
  #dplyr::rename(lon = "longitude"),
  recpointsfill  = "white",
  recpointssize = 3,
  plotsubtitle = expression(paste("Space use by sand tiger sharks ", italic("Carcharias taurus"))),
  pointsincontourssave = paste0(sandtigerSaveDir, "/Scaled/pointsincontours.csv"))  

# sandbar dBBMM --------------------------------------------------------

sandbarSaveDir <- paste0(owd, "/sandbarAllRasters") #new folders for each species and 

### could potentially make this a for loop by making a list of each deployment (or year) and each species
### then it would filter through all combos, creating separate rasters for each 
sandbar_dBBMM <- orstedPosMerged2 %>% 
  filter(common_name_e == "Sandbar") %>% 
  mutate(est = as.POSIXct(est), 
         time = as.POSIXct(time),
         name = as.factor(name)) %>% 
  arrange(est) %>% 
  arrange(name)
str(sandbar_dBBMM)

### When used together, the order of functions would be: movegroup, scaleraster, alignraster if required, plotraster.
### movegroup is nice because it removes individuals that do not have enough positions/detections
movegroup(
  data = sandbar_dBBMM,
  ID = "name",
  Datetime = "est",
  Lat = "latitude",
  Lon = "longitude",
  savedir = sandbarSaveDir) #works with sandbar, sets temporary working directory for individual lever UDs

# Having run the movegroup function example:
scaleraster(path = sandbarSaveDir) #owd = output working directory established early on

# Weighted by number of positions per ID, fewer locations = lower Weighting value = higher final 
# UD values after dividing by Weighting. This scales all IDs up to match the group max.
# I think this is only necessary if there is imbalance in receiver locations but since I'm working with positions not really necessary

#scaleraster(path = tempdir(), weighting = Weighting)

#alginraster combines region-specific group-level UD into a single raster

stationLatsD2 <- stations %>% 
  filter(dep_period == "D1 target") %>% 
  select(latitude) %>% 
  dplyr::rename(lat = "latitude")

stationLonD2 <- stations %>% 
  filter(dep_period == "D1 target") %>% 
  select(longitude) %>% 
  dplyr::rename(lon = "longitude")

plotraster(
  x = paste0(sandbarSaveDir, "/Scaled/All_Rasters_Scaled_Weighted_UDScaled.asc"),
  mapsource = "google",
  maptype = "satellite",
  savedir = paste0(sandbarSaveDir, "/Plot"),
  xlatlon = paste0(sandbarSaveDir, "/Scaled/All_Rasters_Scaled_Weighted_LatLon.asc"),
  locationpoints = sandbar_dBBMM |> dplyr::rename(lat = "latitude", lon = "longitude"),
  #receiverlats = subset(stations, stations$dep_period == "D2 target") |> 
  #dplyr::rename(lat = "latitude"),
  receiverlats = stationLatsD2$lat,
  receiverlons = stationLonD2$lon, 
  #receiverlons = subset(stations, stations$dep_period == "D2 target") |> 
  #dplyr::rename(lon = "longitude"),
  recpointsfill  = "white",
  recpointssize = 3,
  plotsubtitle = expression(paste("Space use by sandbar sharks ", italic("Carcharhinus plumbeus"))),
  pointsincontourssave = paste0(sandbarSaveDir, "/Scaled/pointsincontours.csv"))  

# smooth dogfish dBBMM ----------------------------------------------------

  ### not enough positions right now to run 4/16/2023

