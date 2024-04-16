#Capture Data, Elasmo Size Distribution
#Author: Brianna Cahill
#Purpose: Generate size distribution plots for tagged elasmos
#Input: SunriseWindTaggingData-AcousticTags.csv, data from the SOFOSharks capture data
#Output: Figures summarizing size distribution, csv file with lat/longs of capture locations

# packages and establishing relative pathways -----------------------------
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #Sets wd to where code is
setwd('../Output') #Sets wd to where output file is
owd <- getwd() #Names output file as object owd
setwd('../Data') #Sets wd to where data is
wd <- getwd() #Names data file as object wd

#Packages
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(lubridate)
library(suncalc)
library(dunn.test)
library(FSA)
library(nortest)
library(viridis)
library(tidyverse)
library(stringr)
library(foreach)
library(mgcv)
library(ggpubr)
#library(strptime)
library(car)
library(GGally)
library(ellipse)
library(data.table)
library(stats)
library(MuMIn)
library(ggpubr)
library(lme4)
library(lmerTest)
library(nlme)
library(ggResidpanel)
library(MASS)
library(glmmTMB)
library(emmeans)
library(mgcViz)
library(scales)
library(gratia)
library(ggforce)
library(janitor)
library(gt)
library(glatos)

# formatting data  --------------------------------------------------------

#tagging data
#sunriseTags <- read.csv ("SunriseWindTaggingData - Acoustic Tags.csv", header = TRUE) # i dont think this is necessary
petersonTags <- read.csv ("Sunrise_Orsted_otn_metadata_tagging.xlsx - Tag Metadata.csv", header = TRUE)
petersonTags <- clean_names(petersonTags)

petersonTags$type <- "Animal"
petersonTagsColChange <- petersonTags %>% 
  mutate(transmitter = paste(tag_code_space, tag_id_code, sep = "-")) %>%
  mutate(tag_id = animal_id_floy_tag_id_pit_tag_code_etc) %>%
  mutate(common_name = common_name_e) %>%
  mutate(release_date = utc_release_date_time) %>%
  mutate(tag_expected_life_time_days = est_tag_life) %>%
  mutate(tag_status = "TRUE") %>%
  mutate(measurement = length_m) %>%
  mutate(organization = tag_owner_organization) %>%
  mutate(contact = tag_owner_pi) %>%
  mutate(email = "bradley.peterson@stonybrook.edu")

petersonTagsColChange$common_name[which(petersonTagsColChange$common_name == "Dusky")] <- "Dusky shark"
petersonTagsColChange$common_name[which(petersonTagsColChange$common_name == "Black Sea Bass")] <- "Black sea bass"
petersonTagsColChange$common_name[which(petersonTagsColChange$common_name == "Sand Tiger")] <- "Sand tiger shark"
petersonTagsColChange$common_name[which(petersonTagsColChange$common_name == "Sandbar")] <- "Sandbar shark"
petersonTagsColChange$common_name[which(petersonTagsColChange$common_name == "Smooth Dogfish")] <- "Smooth dogfish"
petersonTagsColChange$common_name[which(petersonTagsColChange$common_name == "Smooth Hammerhead")] <- "Smooth hammerhead shark"
petersonTagsColChange$common_name[which(petersonTagsColChange$common_name == "Winter Skate")] <- "Winter skate"

#adding in the sync tags
syncTags <- read.csv ("Receiever Internal Tag IDs.xlsx - Sheet1.csv", header = TRUE)
syncTags <- clean_names(syncTags)
syncTags$project <- as.factor(syncTags$project)
syncTags$serial_number <- as.factor(syncTags$serial_number)
syncTags$transmit_id <- as.factor(syncTags$transmit_id)
names(syncTags)[names(syncTags) == "transmit_id"] <- "transmitter"
syncTags$type <- "Sync"
syncTags$type <- as.factor(syncTags$type)

#all tags (known and internal sync tags combined for easy filtering/pairing with detections)
knownSyncTags <- dplyr::bind_rows(petersonTagsColChange, syncTags)

#detections
detections2022 <- read.csv ("detection_ORSTED2022_20240106.csv", header = TRUE)
detections2023 <- read.csv ("detection_ORSTED2023_20240106.csv", header = TRUE) #only until early november 2023

detections <- rbind(detections2022, detections2023)
detections <- clean_names(detections)

#detectionsFiltered$detection_timestamp_utc <- as.POSIXct(detectionsFiltered$detection_timestamp_utc, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
#detectionsFiltered$EST <- with_tz(detectionsFiltered$detection_timestamp_utc, "America/New_York")
#detectionsFiltered$receiver_type <- as.factor(detectionsFiltered$receiver_type)
#detectionsFiltered$receiver_sn <- as.factor(detectionsFiltered$receiver_sn)
#detectionsFiltered$Station.Name <- as.factor(detectionsFiltered$Station.Name)
#detectionsFiltered$Transmitter <- as.factor(detectionsFiltered$Transmitter) 
#names(detectionsFiltered)[names(detectionsFiltered) == "Transmitter"] <- "transmitterID"

detections$transmitter <- as.factor(detections$transmitter)
detections$date_and_time_utc <- as.POSIXct(detections$date_and_time_utc, format = "%Y-%m-%d %H:%M:%S", tz = "UTC") 
detections$EST <- with_tz(detections$date_and_time_utc, "America/New_York") #this employs the package lubridate to shift the time
detections$dateEST <- as.Date(detections$EST)

#remove sync tags for ease
mergedDetectionsTags <- merge(detections, unique(knownSyncTags)[, c("tag_id", "transmitter", "tag_sensor_type", "scientific_name", "common_name", 
                                                                    "length_m", "sex", "release_date", "tagger", "type", "contact", "email")], by="transmitter", all.x=TRUE)
summary(as.factor(mergedDetectionsTags$type))
mergedDetectionsTags <- subset(mergedDetectionsTags, type == "Animal")
mergedDetectionsTags$station_name <- as.factor(mergedDetectionsTags$station_name)
nrow(mergedDetectionsTags)

#removing false detections 
mergedDetectionsTags <- mergedDetectionsTags %>% 
  mutate(Transmitter2 = transmitter) %>%
  separate(Transmitter2, c('a', 'b', "transmitter_id"), "-") %>%
  mutate(transmitter_codespace = paste(a, b, sep = "-")) %>%
  separate(receiver, c("receiver_type", "receiver_sn"), "-") %>%
  mutate(detection_timestamp_utc = date_and_time_utc)

mergedDetectionsTags <- false_detections(det=mergedDetectionsTags,
                         tf=3600, #3600 s time threshold
                         show_plot=TRUE)
summary(as.factor(mergedDetectionsTags$passed_filter))
mergedDetectionsTags <- mergedDetectionsTags %>%
  filter(passed_filter == 1)

# tagging humble brags ----------------------------------------------------

test <- petersonTagsColChange %>% 
  filter(tag_sensor_type == "Temperature") %>% 
  group_by(tagger) %>% 
  summarise(count = n())

test2 <- mergedDetectionsTags %>%
  filter(tag_sensor_type == "Temperature") %>% 
  group_by(tagger) %>%
  summarise(animals = n_distinct(tag_id))


# unique individuals ------------------------------------------------------

petersonAnimalsDetected <- mergedDetectionsTags %>% 
  group_by(common_name) %>% 
  summarise(indiv = n_distinct(tag_id))

allAnimalsDetected <- mergedDetectionsTags %>% 
  summarise(indiv = n_distinct(transmitter))

# map stuff ---------------------------------------------------------------

#----- exporting CSV for shark capture locations -----#
tagLocations <- petersonTagsColChange %>%
  filter(tag_sensor_type == "Temperature") #winter skate only has temperature sensor, not depth
write.csv(tagLocations,"tagCaptureLocations.csv", na = "") #replace all NAs with empty cells

# figures ----------------------------------------------------------

#----- size distributions of tagged animals -----#
sizeDist <- petersonTagsColChange %>%
  filter(tag_sensor_type == "Temperature") %>% 
  group_by(common_name, sex) %>%
  summarise(count = n(), 
            mean = mean(length_m), 
            SD = sd(length_m), 
            se = SD/(sqrt(count))) 

sexLabels <- c("F" = "Female", "M" = "Male") #relabeling the facet labels

sizeDistViolinPlot <- test %>%
  filter(tag_sensor_type == "Temperature") %>% 
  filter(common_name != "Smooth hammerhead shark") %>% 
  ggplot(aes(x = common_name, y = length_m, fill = common_name)) + 
  geom_violin(width = 2) +
  geom_sina() + #adds th jitter but confines it to the violin plot
  #geom_text(aes(label = n), data = sizeDist) +
  #facet_wrap(~sex) +
  facet_wrap(~sex, labeller = as_labeller(sexLabels)) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a")) + #dark light dark, using hexadecimal values based on the RdBu palette 
  theme_bw() +
  labs(x= "Common name", y ="Total length (m)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5,), 
        axis.title.x = element_text(size= 18), 
        axis.text.x=element_text(size=16, color="black", angle = 30, hjust = 1),
        axis.title.y = element_text(size= 18),
        axis.text.y=element_text(size=16, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"), 
        legend.position = "none",
        strip.text.x = element_text(size = 16, face = "bold"))
sizeDistViolinPlot
ggsave(paste0(owd,"/","SharkSizeDistributions.png"), width = 19, height = 10)  

#----- abacus plot for tagged animals near orsted receivers (n = 4, date = 5/7/23)-----#
lims <- as.POSIXct(strptime(c("2022-05-15 00:00:00", "2023-11-30 11:59:59"), format = "%Y-%m-%d %H:%M:%S"))

#orstedAnimalDetections <- merge(detections, animalInfo, by = "Transmitter", all = TRUE) 

orstedAnimalDetections <- mergedDetectionsTags %>%
  filter(!is.na(tag_id)) %>%
  filter(!is.na(EST))

check <- mergedDetectionsTags %>%
  filter(tag_sensor_type == "Temperature") %>% 
  group_by(tag_id, common_name, release_date) %>%
  summarise(count = n())

abacusOrsted <- orstedAnimalDetections %>%
  filter(tag_sensor_type == "Temperature") %>% 
  filter(common_name != "Smooth hammerhead shark") %>% 
  ggplot() +
  geom_point(aes(x = EST, y = tag_id, col = common_name), size = 4) +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a")) + 
  geom_point(aes(x = as.POSIXct(release_date, format = "%Y-%m-%d"), y = tag_id), size = 4, shape = 4) +
  geom_vline(aes(xintercept = as.POSIXct("2022-07-30 00:00:00")), linetype ="solid", color = "black", linewidth = 1) +
  facet_grid(common_name ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.placement = "outside") +
  theme_bw() +
  scale_x_datetime(limits=lims, breaks = date_breaks("2 months"), labels=date_format("%b\n%Y")) +
  labs(title= "", x ="Date", y ="Animal Tag ID") + 
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5,),
        axis.title.x = element_text(size=18, face = "bold"), 
        axis.text.x=element_text(size=16, color="black", hjust = 0.5),
        axis.title.y = element_text(size= 18, face = "bold"),
        axis.text.y=element_text(size=8, color="black"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 10, face = "bold"), 
        legend.position = "none")
abacusOrsted
ggsave(paste0(owd,"/","SharkPresenceOrsted.png"), width = 19, height = 10)    
   
#brittney's depth use plot used stat_half eye (for the distribution), geom_points and geom box plot!

#----- abacus plot for tagged animals near orsted receivers (n = 4, date = 5/7/23)-----#
depthUsePlot <- test %>%
  filter(common_name != "Smooth hammerhead shark") %>% 
  filter(tag_sensor_type == "Depth") %>% 
  ggplot(aes(x = common_name, y = sensor_value, fill = common_name)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA) +
  geom_point(
    size = 1.5,
    alpha = .1,
    position = position_jitter(
      seed = 1, width = .1)) + 
  scale_y_reverse() +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a")) + #dark light dark, using hexadecimal values based on the RdBu palette 
  coord_cartesian(clip = "off") + 
  theme_bw() +
  labs(x= "Common name", y ="Depth (m)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5,), 
        axis.title.x = element_text(size= 18), 
        axis.text.x=element_text(size=16, color="black"),
        axis.title.y = element_text(size= 18),
        axis.text.y=element_text(size=16, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"), 
        legend.position = "none",
        strip.text.x = element_text(size = 16, face = "bold"))
depthUsePlot

# omitted ------------------------------------------------------------

#Adjusting Structure of sunriseTags
sunriseTags$Field.Number <- as.factor(sunriseTags$Field.Number)
sunriseTags$Date <- as.Date(sunriseTags$Date, "%m/%d/%y")
sunriseTags$Latitude <- str_remove_all(sunriseTags$Latitude, "N") #fixing the latitude BS
sunriseTags[c('Lat_D', 'Lat_M')] <- str_split_fixed(sunriseTags$Latitude, ' ', 2)
sunriseTags$Lat_D <- as.numeric(sunriseTags$Lat_D)
sunriseTags$Lat_M <- as.numeric(sunriseTags$Lat_M)
sunriseTags$Lat_M[is.na(sunriseTags$Lat_M)] <- 0
sunriseTags$Lat_DD <- sunriseTags$Lat_D + (sunriseTags$Lat_M/60)
sunriseTags$Longitude <- str_remove_all(sunriseTags$Longitude, "W") #fixing the longitude BS
sunriseTags[c('Long_D', 'Long_M')] <- str_split_fixed(sunriseTags$Longitude, ' ', 2)
sunriseTags$Long_D <- as.numeric(sunriseTags$Long_D)
sunriseTags$Long_M <- as.numeric(sunriseTags$Long_M)
sunriseTags$Long_M[is.na(sunriseTags$Long_M)] <- 0
sunriseTags$Long_DD <- (-1*(sunriseTags$Long_D + (sunriseTags$Long_M/60)))
sunriseTags$Bait <- as.factor(sunriseTags$Bait)
sunriseTags$Bait.Type <- as.factor(sunriseTags$Bait.Type)
sunriseTags$species_commonname <- as.factor(sunriseTags$species_commonname)
sunriseTags$sex <- as.factor(sunriseTags$sex)
sunriseTags$PCL <- as.numeric(sunriseTags$PCL)
sunriseTags$FL <- as.numeric(sunriseTags$FL)
sunriseTags$TL <- as.numeric(sunriseTags$TL)
sunriseTags$Girth <- as.numeric(sunriseTags$Girth)
sunriseTags$Tag_SN <- as.factor(sunriseTags$Tag_SN)
sunriseTags$Tag_ID <- as.factor(sunriseTags$Tag_ID)    
    
    
