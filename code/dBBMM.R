
# dBBMM testing -----------------------------------------------------------

library(move)

#---- move testing -----#
testData <- orstedPosMerged2 %>% 
  arrange(est) %>% 
  group_by(full_id) %>% 
  filter(full_id %in% sample(levels(full_id),2))
  as.data.frame()

list <- list(unique(orstedPosMerged2$full_id))
stack<-moveStack(list)
myMoveObject <- moveStack(x=testData$longitude, y=testData$latitude, 
              time=as.POSIXct(testData$time, format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
              proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), 
              data=testData, animal=testData$full_id)
str(myMoveObject)
extent(myMoveObject)

plot(myMoveObject, xlab="Longitude", ylab="Latitude", type="l", pch=16, lwd=0.5)
points(myMoveObject, pch=20, cex=0.5)

ggplot(data = myStackDF, aes(x = location.long, y = location.lat, color = trackId)) + 
  geom_path() + geom_point(size = 0.5) + theme_bw() + coord_cartesian()

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


# new test ----------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(move)
library(plotly)
library(rnaturalearth)
library(sf)
library(MetBrewer)
library(tictoc)

acoustic_data <- orstedPosMerged2 %>% 
  arrange(est) %>% 
  filter(factor(type) == "Orsted") %>% 
  group_by(full_id) %>% 
  filter(n_distinct(est) >= 10) %>% 
  as.data.frame()
  
str(acoustic_data)

glimpse(acoustic_data)
summary(acoustic_data)

acoustic_data_split <- acoustic_data %>% 
  split(.$full_id) #splits data frame by indiv ID

# create "move" object
dat.list <- vector("list", length(acoustic_data_split)) #to store dBBMM results
contours <- vector("list", length(acoustic_data_split)) #to store resulting 50% and 95% UD

# Estimate separately by ID
for (i in 1:length(dat.list)) {
  print(paste("ID", acoustic_data_split[[i]]$full_id[1]))  #print current ID
  
  dat.mov <- move(x = acoustic_data_split[[i]]$longitude, 
                  y = acoustic_data_split[[i]]$latitude, 
                  time = acoustic_data_split[[i]]$time, 
                  data = acoustic_data_split[[i]],
                  proj = CRS("+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs"),
                  animal = acoustic_data_split[[i]]$full_id)
  
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










