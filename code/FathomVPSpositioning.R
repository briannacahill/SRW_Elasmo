#Orsted VPS data analysis
#Author: Brianna Cahill
#Purpose: Import animal position data within the SRW array and create figures
#Input: so many CSVs
#Output: hopefully some figures

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
#library(VTrack) #issues downloading
library(mgcv)
library(ggpubr)
#library(strptime)
library(car)
library(GGally)
library(ellipse)
library(data.table)
library(MuMIn)
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
library(readr)
library(readxl)
library(scattermore)
library(ctmm)
library(sf)
library(ggspatial)
library(ggsflabel)

#devtools::install_github("ctmm-initiative/ctmm")
#devtools::install_github("r-spatial/sf")
#devtools::install_github("yutannihilation/ggsflabel")


# TO DO  -------------------------------------------------------------

# importing data ----------------------------------------------------------

orstedD1T1 <- list.files(path = "positioning/Orsted_D1T1_20221223",
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  lapply(\(x) mutate(x, across(SensorUnit, as.numeric))) %>% #some columns were logical while others were numeric, this allows them to be merged
  bind_rows() %>% 
  as.data.frame()
str(orstedD1T1)

orstedD1T2 <- list.files(path = "positioning/Orsted_D1T2_20230304",
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  lapply(\(x) mutate(x, across(SensorUnit, as.numeric))) %>% #some columns were logical while others were numeric, this allows them to be merged
  bind_rows() %>% 
  as.data.frame()
str(orstedD1T2)

orstedAllPos <- rbind(orstedD1T1, orstedD1T2)

orstedAllPos$EST <- with_tz(orstedAllPos$Time, "America/New_York")

#filter out high error estimates HPEm (HPEm > 10; Bohaboy et al., 2022)
orstedAllPos <- orstedAllPos %>% 
  filter(HPEm < 10 | is.na(HPEm))
summary(orstedAllPos$HPEm)

#need to skip first two empty rows, create column headers, then actually read in the CSV
  #doing this because I'm lazy and I don't feel like creating a CSV from something I already have to use in fathom position
  #probably easiest to use the most up to spec sheet containing tag information
myDeviceCols <- as.character(read_excel("FPspec_T1_20221223.xlsx", sheet = "Devices", skip = 2, n_max = 1, col_names = FALSE)) #need to skip first two empty
devices <- read_excel("FPspec_T1_20221223.xlsx", sheet = "Devices", skip = 2, col_names = myDeviceCols)
devices <- as.data.frame(devices[-1,]) #removes first row (duplicate column header)
str(devices)

#also reading in the orsted tag metadata because we'll need species identifications and maybe other things?
orstedTags <- read.csv ("Sunrise_Orsted_otn_metadata_tagging.xlsx - Tag Metadata.csv", header = TRUE) 
#orstedTags <- clean_names(orstedTags)
#orstedTags <- subset(orstedTags, select = -c(transmitter)) %>% 
  #mutate(Name = animal_id_floy_tag_id_pit_tag_code_etc)

orstedTags <- clean_names(orstedTags) %>% 
  mutate(transmitterID = paste(tag_code_space, tag_id_code, sep = "-")) %>%
  mutate(Name = animal_id_floy_tag_id_pit_tag_code_etc, 
         Type = "Orsted") 

syncTags <- read.csv ("Receiever Internal Tag IDs.xlsx - Sheet1.csv", header = TRUE) %>% 
  mutate(transmitterID = Transmit.ID, 
         Type = "Sync") %>% 
  dplyr::select(-Transmit.ID)

knownSyncTags <- dplyr::bind_rows(orstedTags, syncTags) %>% 
  mutate(FullId = transmitterID)

#merge the positions with the individual info
orstedPosMerged <- merge(orstedAllPos, unique(devices)[, c("FullId", "Name", "SensorType")], by="FullId", all.x=TRUE)
summary(as.factor(orstedPosMerged$Name))
orstedPosMerged2 <- merge(orstedPosMerged, unique(knownSyncTags)[, c("FullId", "tag_sensor_type", "common_name_e", "sex", "Type")], by="FullId", all.x=TRUE)
summary(as.factor(orstedPosMerged2$Name))
head(orstedPosMerged2)

#trying to populate common name column with species for the lobsta and horseshoe crabbos
orstedPosMerged2 <- clean_names(orstedPosMerged2) %>% 
  mutate(common_name_e = case_when(str_detect(name, "HSC") ~ "Horseshoe Crab",
                                   str_detect(name, "L") ~ "Lobster",
                                   TRUE ~ as.character(common_name_e))) 
summary(as.factor(orstedPosMerged2$common_name_e))


# summary and csv things --------------------------------------------------

#cornell tags
cceIndiv <- orstedPosMerged2 %>% 
  filter(common_name_e == "Horseshoe Crab" | common_name_e == "Lobster") %>% 
  group_by(name, full_id) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))  #arranges the tibble in descending order by count
cceIndiv
write.csv(cceIndiv, paste0(owd, "/", "ccePositionCounts.csv"))  

cceFathomPositions <- orstedPosMerged2 %>% 
  filter(common_name_e == "Horseshoe Crab" | common_name_e == "Lobster")
cceFathomPositions
write.csv(cceFathomPositions, paste0(owd, "/", "cceFathomPositions.csv"))  

#HPEm values
syncTagPosCounts <- orstedPosMerged2 %>% 
  filter(type == "Sync") %>% 
  group_by(name, full_id) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))

syncTagPositions <- orstedPosMerged2 %>% 
  filter(type == "Sync")
write.csv(syncTagPositions, paste0(owd, "/", "syncTagPositions.csv"))  

msRequest <- orstedPosMerged2 %>% 
  filter(full_id == "A69-1601-60134" | full_id == "A69-1601-60141" | full_id == "A69-1601-60219")
write.csv(msRequest, paste0(owd, "/", "msRequest_20240318.csv")) 

# figures ----------------------------------------------------------

#need to skip first two empty rows, create column headers, then actually read in the CSV
#doing this because I'm lazy and I don't feel like creating a CSV from something I already have to use in fathom position
myStationCols <- as.character(read_excel("FPspec_T1_20221223.xlsx", sheet = "Stations", skip = 2, n_max = 1, col_names = FALSE)) #need to skip first two empty
stations <- read_excel("FPspec_T1_20221223.xlsx", sheet = "Stations", skip = 2, col_names = myStationCols)
stations <- as.data.frame(stations[-1,]) %>% #removes first row (duplicate column header)
  dplyr::select(Name, Latitude, Longitude, Depth) %>% 
  mutate(Name = as.factor(Name), 
         Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude), 
         Depth = as.numeric(Depth))
str(stations)

sf_receivers <- sf::st_as_sf(stations, coords = c("Longitude", "Latitude"), 
                         crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% 
  st_set_crs(32736) %>% 
  sf::st_transform(32618)
mapview::mapview(sf_receivers)

# dusky MCP -------------------------------------------------------------------

dusky1 <- orstedPositionsMerged2 %>% 
  filter(common_name_e == "Dusky")
str(dusky1)

#test.track <- make_track(dusky1, longitude, latitude, time, crs=32618)
test.track <- make_track(dusky1, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

str(test.track)
class(test.track)
test.track
summary(test.track)
get_crs(test.track)
test.track <- transform_coords(test.track, crs_to=32618)

dat.mcp <- hr_mcp(test.track, levels = c(0.5, 0.95))
dat.mcp
head(dat.mcp$data)
plot(dat.mcp, col = c('red','green','blue'))
get_crs(dat.mcp)

sf_duskyPositions <- sf::st_as_sf(dusky1, coords = c("longitude", "latitude"), 
                                  crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% #"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "EPSG:4979"
  sf::st_transform(32618)

#sf_duskyPositions <- sf::st_as_sf(dusky1, coords = c("longitude", "latitude")) %>% 
  #sf::st_set_crs(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  #sf::st_transform(32618)

confColor <- c(alpha("orange", 0.5), alpha("yellow", 0.2))
library(shades)
figure_duskyT1 <- ggplot() +
  geom_sf(data = dat.mcp$mcp, aes(fill = c(alpha("orange", 0.5), alpha("yellow", 0.2))), size = 0.75, linewidth = 0.9) +
  geom_sf(data = sf_duskyPositions, aes(color = "red"), alpha = 0.4, shape = 16) +
  geom_sf(data = sf_receivers, aes(color = "black"), size = 3.5) +#this should be geom_sf but having CRS issues for some reason
  annotation_scale(location = "bl", width_hint = 0.13) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.45, "in"), pad_y = unit(0.01, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c( st_bbox(sf_duskyPositions)[["xmin"]]-0.001,  st_bbox(sf_duskyPositions)[["xmax"]]+0.001), 
  ylim = c( st_bbox(sf_duskyPositions)[["ymin"]]-0.001,  st_bbox(sf_duskyPositions)[["ymax"]]+0.001))  +
  theme_minimal() +
  scale_fill_manual(values = c(alpha("orange", 0.5), alpha("yellow", 0.2)), labels = c("50%", "95%")) + #these effectively make the legend
  scale_color_manual(values = c("black", "red"), labels = c("Receivers", "Animal positions")) + #these effectively make the legend
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white",
                                        color = "white"),
        plot.background = element_rect(fill = "white",
                                       color = "white")) +
  labs(title = expression(paste("Space use by dusky sharks ", italic("Carcharhinus obscurus"))),
       x =NULL, y = NULL, fill = "", color = "") + #, 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size= 14), 
        axis.text.x=element_text(size=12, color="black", hjust = 0.5),
        axis.title.y = element_text(size= 14),
        axis.text.y=element_text(size=12, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) +
  theme(legend.position= c(0.07, 0.95),
        legend.justification= c("left", "top"),
        legend.spacing.y = unit(-0.1, 'cm'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0))
figure_duskyT1
ggsave(paste0(owd,"/","duskysharks_T1.png"))

# sand tiger MCP -------------------------------------------------------------------

sandtiger1 <- orstedPositionsMerged2 %>% 
  filter(common_name_e == "Sand Tiger")

test.track <- make_track(sandtiger1, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

str(test.track)
class(test.track)
test.track
summary(test.track)
get_crs(test.track)
test.track <- transform_coords(test.track, crs_to=32618)

dat.mcp <- hr_mcp(test.track, levels = c(0.5, 0.95))
dat.mcp
head(dat.mcp$data)
plot(dat.mcp, col = c('red','green','blue'))
get_crs(dat.mcp)

sf_sandtigerPositions <- sf::st_as_sf(sandtiger1, coords = c("longitude", "latitude"), 
                                  crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% #"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "EPSG:4979"
  sf::st_transform(32618)

library(shades)
figure_sandtigerT1 <- ggplot() +
  geom_sf(data = dat.mcp$mcp, aes(fill = c(alpha("orange", 0.5), alpha("yellow", 0.2))), size = 0.75, linewidth = 0.9) +
  geom_sf(data = sf_sandtigerPositions, aes(color = "red"), alpha = 0.4, shape = 16) +
  geom_sf(data = sf_receivers, aes(color = "black"), size = 3.5) +#this should be geom_sf but having CRS issues for some reason
  annotation_scale(location = "bl", width_hint = 0.13) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.45, "in"), pad_y = unit(0.01, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c( st_bbox(sf_sandtigerPositions)[["xmin"]]-0.001,  st_bbox(sf_sandtigerPositions)[["xmax"]]+0.001), 
           ylim = c( st_bbox(sf_sandtigerPositions)[["ymin"]]-0.001,  st_bbox(sf_sandtigerPositions)[["ymax"]]+0.001))  +
  scale_fill_manual(values = c(alpha("orange", 0.5), alpha("yellow", 0.2)), labels = c("50%", "95%")) + #these effectively make the legend
  scale_color_manual(values = c("black", "red"), labels = c("Receivers", "Animal positions")) + #these effectively make the legend
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white",
                                        color = "white"),
        plot.background = element_rect(fill = "white",
                                       color = "white")) +
  labs(title = expression(paste("Space use by sand tiger sharks ", italic("Carcharias taurus"))),
       x =NULL, y = NULL, fill = "", color = "") + #, 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size= 14), 
        axis.text.x=element_text(size=12, color="black", hjust = 0.5),
        axis.title.y = element_text(size= 14),
        axis.text.y=element_text(size=12, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) +
  theme(legend.position= c(0.07, 0.95),
        legend.justification= c("left", "top"),
        legend.spacing.y = unit(-0.1, 'cm'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0))
figure_sandtigerT1
ggsave(paste0(owd,"/","sandtigersharks_T1.png"))

# sandbar MCP -------------------------------------------------------------------

sandbar1 <- orstedPositionsMerged2 %>% 
  filter(common_name_e == "Sandbar")

test.track <- make_track(sandbar1, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

str(test.track)
class(test.track)
test.track
summary(test.track)
get_crs(test.track)
test.track <- transform_coords(test.track, crs_to=32618)

dat.mcp <- hr_mcp(test.track, levels = c(0.5, 0.95))
dat.mcp
head(dat.mcp$data)
plot(dat.mcp, col = c('red','green','blue'))
get_crs(dat.mcp)

sf_sandbarPositions <- sf::st_as_sf(sandbar1, coords = c("longitude", "latitude"), 
                                      crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% #"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "EPSG:4979"
  sf::st_transform(32618)

library(shades)
figure_sandbarT1 <- ggplot() +
  geom_sf(data = dat.mcp$mcp, aes(fill = c(alpha("orange", 0.5), alpha("yellow", 0.2))), size = 0.75, linewidth = 0.9) +
  geom_sf(data = sf_sandbarPositions, aes(color = "red"), alpha = 0.4, shape = 16) +
  #geom_point(data = dusky1, aes(longitude, latitude), color = "red", alpha = 0.4, shape = 16) +
  #geom_point(data = stations, aes(Longitude, Latitude), color = "black", size = 3.5) + 
  geom_sf(data = sf_receivers, aes(color = "black"), size = 3.5) +#this should be geom_sf but having CRS issues for some reason
  annotation_scale(location = "bl", width_hint = 0.13) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.45, "in"), pad_y = unit(0.01, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c( st_bbox(sf_sandbarPositions)[["xmin"]]-400,  st_bbox(sf_sandbarPositions)[["xmax"]]+400), 
           ylim = c( st_bbox(sf_sandbarPositions)[["ymin"]]-400,  st_bbox(sf_sandbarPositions)[["ymax"]]+400))  +
  scale_fill_manual(values = c(alpha("orange", 0.5), alpha("yellow", 0.2)), labels = c("50%", "95%")) + #these effectively make the legend
  scale_color_manual(values = c("black", "red"), labels = c("Receivers", "Animal positions")) + #these effectively make the legend
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white",
                                        color = "white"),
        plot.background = element_rect(fill = "white",
                                       color = "white")) +
  labs(title = expression(paste("Space use by sandbar sharks ", italic("Carcharhinus plumbeus"))),
       x =NULL, y = NULL, fill = "", color = "") + #, 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size= 14), 
        axis.text.x=element_text(size=12, color="black", hjust = 0.5),
        axis.title.y = element_text(size= 14),
        axis.text.y=element_text(size=12, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) +
  theme(legend.position= c(0.07, 0.95),
        legend.justification= c("left", "top"),
        legend.spacing.y = unit(-0.1, 'cm'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0))
figure_sandbarT1
ggsave(paste0(owd,"/","sandbarsharks_T1.png"))

# horseshoe crab MCP -------------------------------------------------------------------

hsc1 <- orstedPositionsMerged2 %>% 
  filter(common_name_e == "Horseshoe Crab")

test.track <- make_track(hsc1, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

str(test.track)
class(test.track)
test.track
summary(test.track)
get_crs(test.track)
test.track <- transform_coords(test.track, crs_to=32618)

dat.mcp <- hr_mcp(test.track, levels = c(0.5, 0.95))
dat.mcp
head(dat.mcp$data)
plot(dat.mcp, col = c('red','green','blue'))
get_crs(dat.mcp)

sf_hscPositions <- sf::st_as_sf(hsc1, coords = c("longitude", "latitude"), 
                                      crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% #"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "EPSG:4979"
  sf::st_transform(32618)

library(shades)
figure_hscT1 <- ggplot() +
  geom_sf(data = dat.mcp$mcp, aes(fill = c(alpha("orange", 0.5), alpha("yellow", 0.2))), size = 0.75, linewidth = 0.9) +
  geom_sf(data = sf_hscPositions, aes(color = "red"), alpha = 0.4, shape = 16) +
  #geom_point(data = dusky1, aes(longitude, latitude), color = "red", alpha = 0.4, shape = 16) +
  #geom_point(data = stations, aes(Longitude, Latitude), color = "black", size = 3.5) + 
  geom_sf(data = sf_receivers, aes(color = "black"), size = 3.5) +#this should be geom_sf but having CRS issues for some reason
  annotation_scale(location = "bl", width_hint = 0.13) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.45, "in"), pad_y = unit(0.01, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c( st_bbox(sf_hscPositions)[["xmin"]]-0.001,  st_bbox(sf_hscPositions)[["xmax"]]+0.001), 
           ylim = c( st_bbox(sf_hscPositions)[["ymin"]]-0.001,  st_bbox(sf_hscPositions)[["ymax"]]+0.001))  +
  scale_fill_manual(values = c(alpha("orange", 0.5), alpha("yellow", 0.2)), labels = c("50%", "95%")) + #these effectively make the legend
  scale_color_manual(values = c("black", "red"), labels = c("Receivers", "Animal positions")) + #these effectively make the legend
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white",
                                        color = "white"),
        plot.background = element_rect(fill = "white",
                                       color = "white")) +
  labs(title = expression(paste("Space use by horseshoe crabs ", italic("Limulus polyphemus"))),
       x =NULL, y = NULL, fill = "", color = "") + #, 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size= 14), 
        axis.text.x=element_text(size=12, color="black", hjust = 0.5),
        axis.title.y = element_text(size= 14),
        axis.text.y=element_text(size=12, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) +
  theme(legend.position= c(0.07, 0.95),
        legend.justification= c("left", "top"),
        legend.spacing.y = unit(-0.1, 'cm'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0))
figure_hscT1
ggsave(paste0(owd,"/","horseshoecrabs_T1.png"))

# sync tags MCP -------------------------------------------------------------------

summary(as.factor(orstedPositionsMerged2$type))
syncTags1 <- orstedPositionsMerged2 %>% 
  filter(type == "Sync")

test.track <- make_track(syncTags1, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

str(test.track)
class(test.track)
test.track
summary(test.track)
get_crs(test.track)
test.track <- transform_coords(test.track, crs_to=32618)

dat.mcp <- hr_mcp(test.track, levels = c(0.5, 0.95))
dat.mcp
head(dat.mcp$data)
plot(dat.mcp, col = c('red','green','blue'))
get_crs(dat.mcp)

sf_syncPositions <- sf::st_as_sf(syncTags1, coords = c("longitude", "latitude"), 
                                 crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% #"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "EPSG:4979"
  sf::st_transform(32618)

library(shades)
figure_syncT1 <- ggplot() +
  #geom_sf(data = dat.mcp$mcp, aes(fill = c(alpha("orange", 0.5), alpha("yellow", 0.2))), size = 0.75, linewidth = 0.9) +
  geom_sf(data = sf_syncPositions, aes(color = as.factor(full_id)), alpha = 0.4, shape = 16) +
  #geom_point(data = dusky1, aes(longitude, latitude), color = "red", alpha = 0.4, shape = 16) +
  #geom_point(data = stations, aes(Longitude, Latitude), color = "black", size = 3.5) + 
  geom_sf(data = sf_receivers, aes(color = "black"), size = 3.5) +#this should be geom_sf but having CRS issues for some reason
  annotation_scale(location = "bl", width_hint = 0.13) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.45, "in"), pad_y = unit(0.01, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c( st_bbox(sf_syncPositions)[["xmin"]]-0.001,  st_bbox(sf_syncPositions)[["xmax"]]+0.001), 
           ylim = c( st_bbox(sf_syncPositions)[["ymin"]]-0.001,  st_bbox(sf_syncPositions)[["ymax"]]+0.001))  +
  scale_fill_manual(values = c(alpha("orange", 0.5), alpha("yellow", 0.2)), labels = c("50%", "95%")) + #these effectively make the legend
  #scale_color_manual(values = c("black", "red"), labels = c("Receivers", "Sync tag positions")) + #these effectively make the legend
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white",
                                        color = "white"),
        plot.background = element_rect(fill = "white",
                                       color = "white")) +
  labs(title = "Positions of sync tags within array",
       x =NULL, y = NULL, fill = "", color = "") + #, 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size= 14), 
        axis.text.x=element_text(size=12, color="black", hjust = 0.5),
        axis.title.y = element_text(size= 14),
        axis.text.y=element_text(size=12, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) +
  theme(legend.position= c(0.07, 0.95),
        legend.justification= c("left", "top"),
        legend.spacing.y = unit(-0.1, 'cm'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0))
figure_syncT1
ggsave(paste0(owd,"/","syncTags_T1.png"))


# checking on missing receivers -------------------------------------------

############ working on this in untitled 3 (unless it gets accidentally deleted)
summary(as.factor(orstedPositionsMerged2$Type))
syncTags1 <- orstedPositionsMerged2 %>% 
  filter(Type == "Sync")
summary(as.factor(syncTags1$Name))
summary(syncTags1$HPEm)

missingSyncTags <- orstedPositionsMerged2 %>% 
  filter(Type == "Sync") %>% 
  filter(Name == "R1C1" | Name == "R1C2" | Name == "R1C3" | Name == "R1C4" | Name == "R1C5" | Name == "R1C6" | Name == "R1C8" | Name == "R2C6" |
           Name == "R4C5" | Name == "R4C7") 

lostStationsT1 <- c("R1C1", "R1C2", "R1C3", "R1C4", "R1C5", 
                    "R1C6", "R1C8", "R2C6", "R4C5", "R4C7")  
lostStationsT1.Lat <- c(40.731390, 40.732400, 40.733390, 40.734490, 40.735500, 
                        40.736500, 40.738590, 40.733510, 40.726530, 40.728580)
lostStationsT1.Long <- c(-72.842100, -72.839070, -72.836200, -72.833170, -72.830180, 
                         -72.827290, -72.821270, -72.826290, -72.824220, -72.821290)

lostStationCoords <- cbind(lostStationsT1, lostStationsT1.Lat, lostStationsT1.Long) %>% 
  as.data.frame() %>% 
  mutate(lostStationsT1.Lat = as.numeric(lostStationsT1.Lat), 
         lostStationsT1.Long = as.numeric(lostStationsT1.Long))

test <- ggplot() +
  geom_point(data = missingSyncTags, aes(x = Longitude, y = Latitude, color = Name)) +
  geom_point(data = stations, aes(x = Longitude, y = Latitude), color = "black") +
  geom_point(data = lostStationCoords, aes(x = lostStationsT1.Long, y = lostStationsT1.Lat), color = "black", shape = 4, size = 4) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white",
                                        color = "white"),
        plot.background = element_rect(fill = "white",
                                       color = "white")) +
  labs(title = expression(paste("B)")),
       x =NULL, y = NULL, fill = "", color = "Stations") + #, 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0), 
        axis.title.x = element_text(size= 14), 
        axis.text.x=element_text(size=12, color="black", hjust = 0.5),
        axis.title.y = element_text(size= 14),
        axis.text.y=element_text(size=12, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) +
  theme(legend.position = "none")
    #legend.justification= "none",
    #legend.spacing.y = unit(-0.1, 'cm'))
#legend.margin=margin(0,0,0,0),
#legend.box.margin=margin(0,0,0,0))
test

test2 <- ggplot() +
  geom_point(data = missingSyncTags, aes(x = Time, y = Name, color = Name)) +
  scale_y_discrete(limits=rev) +
  scale_x_datetime(date_breaks = "1 month", date_labels ="%b\n%Y") + #date_minor_breaks = "1 week", 
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white",
                                        color = "white"),
        plot.background = element_rect(fill = "white",
                                       color = "white")) +
  labs(title = expression(paste("A)")),
       x = "Date", y = "Stations", fill = "", color = "Stations") + 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0), 
        axis.title.x = element_text(size= 14), 
        axis.text.x=element_text(size=12, color="black", hjust = 0.5),
        axis.title.y = element_text(size= 14),
        axis.text.y=element_text(size=12, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) +
  theme(legend.position = "none")
test2

missingReceivers <- ggarrange(test2, test,
                                  ncol = 2, nrow = 1)
missingReceivers
#annotate_figure(allSpeciesDetections, left = text_grob("Total detections", rot = 90, size = 20, face = "bold"),
                #bottom = text_grob("Month", size = 20, face = "bold"))
ggsave(paste0(owd,"/","missingReceiverPositions.png"), width = 19, height = 10) 

# notes -------------------------------------------------------------------

#COA in VTrack
#igraph package for network analyses
  #approach links the movements of each tagged individual to each visited area (Dakity or Manglar Bay) by an edge (arrow) and is weighted by the 
    #number of movements/detections, thus providing an indicator of the use of each region
  #produce bipartite graphs
  #Brownscombe et al., () paper did this with bonefish comparing movements between bays, could use this for looking at movement between ACOE borrows
  #state-space models using aniMotum R package (Jonsen et al., 2023)







