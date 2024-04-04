
# dBBMM testing -----------------------------------------------------------

library(move)
library(actel) #not currently using
library(RSP) #not currently using
library(movegroup) #not currently using

#---- move testing -----#
testData <- orstedPosMerged2 %>% 
  arrange(est) %>% 
  group_by(full_id) %>% 
  as.data.frame()
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

# Create a moveGroup object
movegroup <- movegroup(acoustic_data, lon = "longitude", lat = "latitude", time = "est", id = "full_id")

# Prepare data for RSP
rsp_data <- move2dt(movegroup)

# Estimate utilization distribution (UD)
UD <- UtilizationDistribution(rsp_data)

# Fit dBBMM
fit <- FitDBBMM(rsp_data, UD)

# Summary of the dBBMM results
summary(fit)

# Plot dBBMM results
plot(fit)
