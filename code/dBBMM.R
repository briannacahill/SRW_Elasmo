
# dBBMM testing -----------------------------------------------------------

library(move)

#---- move testing -----#
duskyD1 <- orstedPosMerged2 %>% 
  filter(common_name_e == "Dusky") %>% 
  filter(dep_period == "Deployment 1") %>% 
  arrange(est) %>% 
  arrange(name) %>% 
  as.data.frame()
str(duskyD1)

test <- duskyD1 %>% 
  group_by(name) %>% 
  summarise(count = n())
test # GM2211 has the most positions


duskyD1 <- orstedPosMerged2 %>% 
  filter(common_name_e == "Dusky") %>% 
  filter(dep_period == "Deployment 1") %>% 
  mutate(est = as.POSIXct(est), 
         time = as.POSIXct(time),
         name = as.factor(name)) %>% 
  arrange(est) %>% 
  arrange(name) %>% 
  filter(n_distinct(est) >= 25) %>% 
  as.data.frame()
str(duskyD1)

#test.track <- make_track(dusky1, longitude, latitude, time, crs=32618)
test.track <- make_track(duskyD1, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

dfGM2211move <-  move(x=duskyD1$longitude, y=duskyD1$latitude, 
                      time=duskyD1$time, 
                      proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), 
                      data=duskyD1, animal=duskyD1$name, sensor=duskyD1$sensor_value)

str(dfGM2211move)
extent(dfGM2211move)
namesIndiv(dfGM2211move)

plot(dfGM2211move, xlab="Longitude", ylab="Latitude", type="l", pch=16, lwd=0.5)
points(dfGM2211move, pch=20, cex=0.5)

myStackDF <- as.data.frame(dfGM2211move)
ggplot(data = myStackDF, aes(x = longitude, y = latitude, color = trackId)) + 
  geom_path() + geom_point(size = 0.5) + theme_bw() + coord_cartesian()

dfGM2211moveCorr <- corridor(dfGM2211move)
plot(dfGM2211moveCorr, type="l", xlab="Longitude", ylab="Latitude", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

dfGM2211move.prj <- spTransform(dfGM2211move, center=TRUE)
min(timeLag(x = dfGM2211move,units = "mins")) 

dBB.dfGM2211 <- brownian.bridge.dyn(dfGM2211move.prj, ext= 4, raster=100, location.error=20) #ext default is .3 (30% expansion on raster area)
$bbox[1,]
plot(dBB.dfGM2211, main="dBBMM")

## calculate the dynamic brownian motion variance of the gappy track
dbbv <- brownian.motion.variance.dyn(dfGM2211move.prj, location.error=20, window.size=31, margin=11)

## the intended GPS fix rate of leroy was 15min, so we will ignore for example all segments that have a larger time lag than 5hours. The 'dBMvariance' object resulting from the function above, contains the slot '@interest' in which those segments marked as FALSE won't be included in the calculation of the dBBMM. Therefore we set all segments with time lag larger than 300mins to false
dbbv@interest[timeLag(dfGM2211move.prj,"mins")>300] <- FALSE


#----- switching over to leroy because it follows the tutorial -----#
leroyRed <- leroy[1:200] # reducing dataset for example
leroy.prj <- spTransform(leroyRed, center=TRUE) # center=T: the center of the coordinate system is the center of the track. Units are in meters
summary(timeLag(leroyRed,"mins")) 

dBB.leroy <- brownian.bridge.dyn(leroy.prj, ext=.85, raster=100, location.error=20)
plot(dBB.leroy, main="dBBMM")

UDleroy <- getVolumeUD(dBB.leroy)

par(mfrow=c(1,2))
plot(UDleroy, main="UD")

## also a contour can be added
plot(UDleroy, main="UD and contour lines")
contour(UDleroy, levels=c(0.5, 0.95), add=TRUE, lwd=c(0.5, 0.5), lty=c(2,1))

par(mfrow=c(1,3))

## mantaining the lower probabilities
ud95 <- UDleroy
ud95[ud95>.95] <- NA
plot(ud95, main="UD95")

## or extracting the area with a given probability, where cells that belong to the given probability will get the value 1 while the others get 0
ud95 <- UDleroy<=.95
plot(ud95, main="UD95")

ud50 <- UDleroy<=.5
plot(ud50, main="UD50")

## creating a gappy data set
leroyWithGap <- leroy[-c(50:500,550:850)]
leroyWithGap_p <- spTransform(leroyWithGap, center=TRUE)

## calculate the dBBMM with the default extent gives an error that it is too small
dbb <- brownian.bridge.dyn(leroyWithGap_p, raster=100, location.error=20)

## making the extent bigger seems to solve the problem
dbb <- brownian.bridge.dyn(leroyWithGap_p, raster=100, location.error=20, ext=4)

## but than the UD is not very informative
ud <- getVolumeUD(dbb)

par(mfrow=c(1,2))
plot(ud, main="UD")
contour(ud, levels=c(0.5, 0.95), add=TRUE, lwd=c(0.5, 0.5), lty=c(2,1))

plot(ud, main="UD with locations")
points(leroyWithGap_p, col="red",  cex=.5, pch=20)
contour(ud, levels=c(0.5, 0.95), add=TRUE, lwd=c(0.5, 0.5), lty=c(2,1))

## calculate the dynamic brownian motion variance of the gappy track
dbbv <- brownian.motion.variance.dyn(leroyWithGap_p, location.error=20, window.size=31, margin=11)

## the intended GPS fix rate of leroy was 15min, so we will ignore for example all segments that have a larger time lag than 5hours. The 'dBMvariance' object resulting from the function above, contains the slot '@interest' in which those segments marked as FALSE won't be included in the calculation of the dBBMM. Therefore we set all segments with time lag larger than 300mins to false
dbbv@interest[timeLag(leroyWithGap_p,"mins")>300] <- FALSE
  # cant get this part to work I think because all indiv are included in this

## then we use the 'dBMvariance' object to calculate the dBBMM
dbb.corrected <- brownian.bridge.dyn(dbbv, raster=100, ext=.45,location.error=20)

## now the UD makes more sense
ud.corrected <- getVolumeUD(dbb.corrected)

par(mfrow=c(1,2))
plot(ud.corrected, main="UD")
contour(ud.corrected, levels=c(0.5, 0.95), add=TRUE, lwd=c(0.5, 0.5), lty=c(2,1))

plot(ud.corrected, main="UD with locations")
points(leroyWithGap_p, col="red", cex=.5, pch=20)
contour(ud.corrected, levels=c(0.5, 0.95), add=TRUE, lwd=c(0.5, 0.5), lty=c(2,1))



########################## got up until here to work ##########################

library(movegroup)

data(TracksCleaned)
View(TracksCleaned)

# 1. movegroup
movegroup(data = TracksCleaned, ID = "Shark", Datetime = "Datetime", Lat = "Lat", Lon = "Lon", savedir = "~/Downloads")
str(test)

# 2. scaleraster

# 3. alignraster (if required)

# 4. plotraster

########################## stuff that I haven't really gotten to work ##########################
tic()
dat.list <- move::brownian.bridge.dyn(object = dfGM2211move, raster = 0.25, location.error=12,
                                      margin = 9, window.size = 29, ext = rast.ext)
toc()


## Extract 50 and 95% contours of space-use
res <- raster2contour(dat.list[[i]], levels = c(0.5, 0.95))

contours[[i]] <- st_as_sf(res)

if (st_geometry_type(contours[[i]])[1] == 'LINESTRING') {
  contours[[i]] <- st_cast(contours[[i]], "POLYGON")
} else {
  contours[[i]] <- st_cast(contours[[i]], "MULTIPOLYGON") %>%
    st_make_valid()  #fixes issue w/ negative areas being calculated
}






data2 <- spTransform(dfGM2211move, CRSobj="+proj=aeqd +ellps=WGS84", center=TRUE)

dbbmm <- brownian.bridge.dyn(object=data2, location.error=201, dimSize=125, ext=1.2, 
                             time.step=2, margin=31) #, location.error=12, dimSize=125, ext=1.2, time.step=2, margin=15

plot(dbbmm)

#ggplot(data = dfGM2211move, aes(x = location.long, y = location.lat, color = trackId)) + 
#geom_path() + geom_point(size = 0.5) + theme_bw() + coord_cartesian()

library(crawl)

dbbmm_fit <- crawl::crawl(dfGM2211move)


#----- movegroup, rsp and actel testing -----#

# Load required libraries
library(movegroup)
library(actel)
library(RSP)

# Sample data - replace this with your actual data
# Assuming you have a data frame called 'acoustic_data' with columns 'timestamp', 'receiver_id', 'animal_id', 'latitude', and 'longitude'
# 'timestamp' is the time of each detection
# 'receiver_id' is the identifier for the acoustic receiver
# 'animal_id' is the identifier for the tagged animal
# 'latitude' and 'longitude' are the coordinates
# Example data:
# acoustic_data <- data.frame(timestamp = c("2024-01-01 12:00:00", "2024-01-01 12:10:00", "2024-01-01 12:20:00"),
#                             receiver_id = c(1, 1, 2),
#                             animal_id = c("A", "B", "A"),
#                             latitude = c(40.0, 40.1, 40.2),
#                             longitude = c(-70.0, -70.1, -70.2))

acoustic_data <- orstedPosMerged2 %>% 
  arrange(est) %>% 
  group_by(full_id) %>% 
  as.data.frame()
str(acoustic_data)

test <- orstedPosMerged2 %>%
  group_by(common_name_e) %>% 
  summarise(count = n(), 
            indiv = length(unique(full_id)))
test

sandtiger_data <- orstedPosMerged2 %>% 
  filter(common_name_e == "Sand Tiger") %>% 
  arrange(est) %>% 
  group_by(full_id) %>% 
  as.data.frame()
str(sandtiger_data)

# Create a moveGroup object
movegroup <- movegroup(sandbar_data, lon = "longitude", lat = "latitude", time = "est", id = "full_id")
movegroup <- movegroup(sandbar_data, lon = longitude, lat = latitude, time = est, id = full_id)

#When used together, the order of functions would be: movegroup, scaleraster, alignraster if required, plotraster.
test <- movegroup(
  data = sandtiger_data,
  ID = "full_id",
  Datetime = "est",
  Lat = "latitude",
  Lon = "longitude",
  savedir = tempdir()) #works with sandbar, sets temporary working directory for individual lever UDs
test

# Having run the movegroup function example:
scaleraster(path = tempdir())

# Weighted by number of positions per ID, fewer locations = lower Weighting value = higher final 
# UD values after dividing by Weighting. This scales all IDs up to match the group max.
Weighting <- sandtiger_data |>
  dplyr::group_by(full_id) |>
  dplyr::summarise(N = n()) |> 
  dplyr::filter(N > 23) |> 
  dplyr::mutate(N = N / max(N, na.rm = TRUE)) |> 
  dplyr::pull(N)

scaleraster(path = tempdir(), weighting = Weighting)


# Josh Cullen Youtube video testing, for loop addressing indiv  ----------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(move)
library(plotly)
library(rnaturalearth)
library(sf)
library(MetBrewer)
library(tictoc)

acoustic_data <- orstedPosMerged2 %>% 
  filter(common_name_e == "Dusky") %>% 
  filter(dep_period == "Deployment 1") %>% 
  mutate(est = as.POSIXct(est), 
         time = as.POSIXct(time),
         name = as.factor(name)) %>% 
  arrange(est) %>% 
  arrange(name) %>% 
  filter(n_distinct(est) >= 10) %>% 
  as.data.frame()

  

str(acoustic_data)

glimpse(acoustic_data)
summary(acoustic_data)

#acoustic_data_split <- acoustic_data %>% 
  #split(.$name) #splits data frame by indiv ID

# create "move" object
dat.list <- vector("list", length(acoustic_data)) #to store dBBMM results
contours <- vector("list", length(acoustic_data)) #to store resulting 50% and 95% UD

# Estimate separately by ID
for (i in 1:length(dat.list)) {
  print(paste("ID", acoustic_data[[i]]$name[1]))  #print current ID
  
  dat.mov <- move(x = acoustic_data[[i]]$longitude, 
                  y = acoustic_data[[i]]$latitude, 
                  time = acoustic_data[[i]]$time, 
                  proj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                  data = acoustic_data[[i]],
                  animal = acoustic_data[[i]]$name)
  
  # Conditionally define extent; necessary for turtles that don't migrate
  x.ext <- diff(dat.mov$bbox[1,])
  rast.ext <- ifelse(x.ext < 100, 4, 0.3)
  
  
  ## Run dBBMM (this will take a little while to run)
  
  tic()
  dat.list[[i]] <- brownian.bridge.dyn(object = dat.mov, raster = 100, location.error=20, ext = rast.ext)
  toc()
  
  ## Extract 50 and 95% contours of space-use
  res <- raster2contour(dat.list[[i]], levels = c(0.5, 0.95))
  
  contours[[i]] <- st_as_sf(res)
  
  if (st_geometry_type(contours[[i]])[1] == 'LINESTRING') {
    contours[[i]] <- st_cast(contours[[i]], "POLYGON")
  } else {
    contours[[i]] <- st_cast(contours[[i]], "MULTIPOLYGON") %>%
      st_make_valid()  #fixes issue w/ negative areas being calculated
  }
  
}

# Estimate separately by ID
for (i in 1:length(dat.list)) {
  print(paste("ID", acoustic_data_split[[i]]$name[1]))  #print current ID
  
  dat.mov <- move(x = acoustic_data_split[[i]]$longitude, 
                  y = acoustic_data_split[[i]]$latitude, 
                  time = acoustic_data_split[[i]]$time, 
                  data = acoustic_data_split[[i]],
                  proj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                  animal = acoustic_data_split[[i]]$name)
  
  # Conditionally define extent; necessary for turtles that don't migrate
  x.ext <- diff(dat.mov$bbox[1,])
  rast.ext <- ifelse(x.ext < 100, 3, 0.3)
  
  
  ## Run dBBMM (this will take a little while to run)
  
  tic()
  dat.list[[i]] <- brownian.bridge.dyn(object = dat.mov, raster = 0.25, location.error=12,
                                       margin = 9, window.size = 29, ext = rast.ext)
  toc()
  
  
  ## Extract 50 and 95% contours of space-use
  res <- raster2contour(dat.list[[i]], levels = c(0.5, 0.95))
  
  contours[[i]] <- st_as_sf(res)
  
  if (st_geometry_type(contours[[i]])[1] == 'LINESTRING') {
    contours[[i]] <- st_cast(contours[[i]], "POLYGON")
  } else {
    contours[[i]] <- st_cast(contours[[i]], "MULTIPOLYGON") %>%
      st_make_valid()  #fixes issue w/ negative areas being calculated
  }
  
}

# subsetting the data

acoustic_dataSUBSET <- orstedPosMerged2 %>% 
  arrange(est) %>% 
  filter(full_id == "A69-9004-11528") %>% 
  as.data.frame()
str(acoustic_dataSUBSET)

dat.mov <- move(x = acoustic_data_split[[i]]$longitude, 
                y = acoustic_data_split[[i]]$latitude, 
                time = acoustic_data_split[[i]]$time, 
                data = acoustic_data_split[[i]],
                proj = CRS("+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs"),
                animal = acoustic_data_split[[i]]$full_id)

tic()
test <- brownian.bridge.dyn(object = dat.mov, raster = 0.25,  location.error=12, dimSize=125, ext=1.2, 
                            time.step=2, margin=15)
toc()

data(leroy)

# RSP heat map ---------------------------------------------------------
library(remotes)
library(RSP)
browseVignettes("RSP")
#install_github("YuriNiella/RSP", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

States_shapefile-shp

land <- actel::shapeToRaster(shape = wd, "/States_shapefile-shp/States_shapefile.shp"), 
                              size = 0.0001, buffer = 0.05) 

# Create a transition layer with 8 directions
tl <- actel::transitionLayer(x = water, directions = 8)

# Import example output from actel::explore() 
data(input.example) 

# Run RSP analysis
rsp.data <- runRSP(input = input.example, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude")










