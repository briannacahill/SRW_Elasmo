#Orsted VPS data analysis
#Author: Brianna Cahill
#Purpose: Import consolidated animal and sync tag positions and then generate figures
#Input: only a few CSVs
#Output: hopefully some figures

# packages and establishing relative pathways -----------------------------

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
library(amt)
library(shades)

#devtools::install_github("ctmm-initiative/ctmm")
#devtools::install_github("r-spatial/sf")
#devtools::install_github("yutannihilation/ggsflabel")

# relative pathways
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #Sets wd to where code is
setwd('../Output') #Sets wd to where output file is
owd <- getwd() #Names output file as object owd
setwd('../Data') #Sets wd to where data is
wd <- getwd() #Names data file as object wd

# reading in CSVs ---------------------------------------------------------

sbuTags <- read.csv("sbuFathomPositions.csv", header = TRUE) 
syncTags <- read.csv("syncTagPositions.csv", header = TRUE) 

orstedPosMerged2 <- rbind(sbuTags, syncTags)

# error comparison -----------------------------------------------------

summary(orstedPosMerged2$hp_em)
summary(orstedPosMerged2$hp_es)
summary(orstedPosMerged2$rmse)
summary(as.factor(orstedPosMerged2$dep_period))

test <- orstedPosMerged2 %>% 
  #filter(hp_em < 10) %>% 
  group_by(dep_period) %>% 
  summarise(count = n(), 
            as_tibble_row(quantile(hp_es, na.rm = TRUE)), 
            mean = mean(hp_es, na.rm = TRUE), 
            sd = sd(hp_es, na.rm = TRUE))
test

kruskal.test(hp_em ~ dep_period, data = orstedPosMerged2)
kruskal.test(hp_em ~ dep_period, data = orstedPosMerged2[orstedPosMerged2$hp_em < 10,])

# i think this is showing distance, isn't the HPE value that we want
plotHPEm <- orstedPosMerged2 %>% 
  filter(hp_em < 10) %>% 
  group_by(dep_period) %>%
  ggplot() +
  geom_violin(aes(x = dep_period, y = hp_em)) +
  theme_minimal() +
  labs(title = NULL, x =NULL, y = "HPEm")

plotHPEs <- orstedPosMerged2 %>% 
  filter(hp_es < 50) %>% 
  group_by(dep_period) %>%
  ggplot() +
  geom_violin(aes(x = dep_period, y = hp_es)) +
  theme_minimal() +
  labs(title = NULL, x =NULL, y = "HPE") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size= 14), 
        axis.text.x=element_text(size=12, color="black", hjust = 0.5),
        axis.title.y = element_text(size= 14),
        axis.text.y=element_text(size=12, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"))

plotRMSE <- orstedPosMerged2 %>% 
  filter(rmse < 10) %>% 
  group_by(dep_period) %>%
  ggplot() +
  geom_violin(aes(x = dep_period, y = rmse)) +
  theme_minimal() +
  labs(title = NULL, x =NULL, y = "RMSE") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size= 14), 
        axis.text.x=element_text(size=12, color="black", hjust = 0.5),
        axis.title.y = element_text(size= 14),
        axis.text.y=element_text(size=12, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"))

errorComp <- ggarrange(plotHPEs, plotRMSE, 
                       ncol = 2, nrow = 1)
errorComp
ggsave(paste0(owd,"/","errorComp_deployment.png"))

# figures ----------------------------------------------------------


# station info -------------------------------------------------

stations <- read.csv("sunrise_coords.csv", header = TRUE) %>% 
  mutate(dep_period = as.factor(dep_period), 
         name = as.factor(name),
         status = as.factor(status))

sf_stations <- sf::st_as_sf(stations, coords = c("longitude", "latitude"), 
                             crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% 
  #st_set_crs(32736) %>% 
  sf::st_transform(32618)

# deployment time periods -------------------------------------------------

D1start <- as.POSIXct("2022-07-30 12:00:00")
D1end <- as.POSIXct("2023-05-17 08:00:00")
D2start <- as.POSIXct("2023-05-17 15:55:00")
D2end <- as.POSIXct("2023-10-27 08:00:00")

# dusky MCP -------------------------------------------------------------------

dusky1 <- orstedPosMerged2 %>% 
  filter(common_name_e == "Dusky") %>% 
  mutate(est = as.POSIXct(est), 
         time = as.POSIXct(time),
         name = as.factor(name), 
         year = format(as.POSIXct(est), format = "%Y")) %>% 
  arrange(est) %>% 
  arrange(name)

du.track <- make_track(dusky1, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

str(du.track)
class(du.track)
du.track
summary(du.track)
get_crs(du.track)
du.track <- transform_coords(du.track, crs_to=32618)

du.mcp <- hr_mcp(du.track, levels = c(0.5, 0.95))
du.mcp
head(du.mcp$data)
plot(du.mcp, col = c('red','green','blue'))
get_crs(du.mcp)

sf_duskyPositions <- sf::st_as_sf(dusky1, coords = c("longitude", "latitude"), 
                                    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% #"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "EPSG:4979"
  sf::st_transform(32618)


figure_duskyT1 <- #ggmap(SRWbasemap) +
  ggplot() +
  #minimum convex polygons
  geom_sf(data = du.mcp$mcp, 
          mapping = aes(fill = c(alpha("orange", 0.5), alpha("yellow", 0.2))), 
          size = 0.75, 
          linewidth = 0.9, 
          inherit.aes = FALSE) +
  #geom_sf(data = subset(sf_sandbarPositions, sf_sandbarPositions$dep_period), 
  #aes(color = "red"), alpha = 0.4, shape = 16) +
  #sandbar positions colored by year
  geom_sf(data = sf_duskyPositions, 
          mapping = aes(color = year), 
          alpha = 0.6, 
          shape = 16, 
          inherit.aes = FALSE) +
  #receiver stations showing target locations
  geom_sf(data = subset(sf_stations, sf_receivers$dep_period == "D1 target"), 
          mapping = aes(color = "black"), 
          size = 3.5, 
          inherit.aes = FALSE) +
  coord_sf(crs = st_crs(32618)) + #3857
  annotation_scale(location = "bl", width_hint = 0.13) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.45, "in"), pad_y = unit(0.01, "in"),
                         style = north_arrow_fancy_orienteering) +
  #coord_sf(xlim = c( st_bbox(sf_sandbarPositions)[["xmin"]]-400,  st_bbox(sf_sandbarPositions)[["xmax"]]+400), 
  #ylim = c( st_bbox(sf_sandbarPositions)[["ymin"]]-400,  st_bbox(sf_sandbarPositions)[["ymax"]]+400))  +
  scale_fill_manual(values = c(alpha("orange", 0.5), alpha("yellow", 0.2)), labels = c("50%", "95%")) + #these effectively make the legend
  scale_color_manual(values = c("red", "blue", "black"), 
                     labels = c("2022 positions", "2023 positions", "Target receiver locations")) + #these effectively make the legend
  theme_minimal() +
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

sandtiger1 <- orstedPosMerged2 %>% 
  filter(common_name_e == "Sand Tiger") %>% 
  mutate(est = as.POSIXct(est), 
         time = as.POSIXct(time),
         name = as.factor(name), 
         year = format(as.POSIXct(est), format = "%Y")) %>% 
  arrange(est) %>% 
  arrange(name)

st.track <- make_track(sandtiger1, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

str(st.track)
class(st.track)
st.track
summary(st.track)
get_crs(st.track)
st.track <- transform_coords(st.track, crs_to=32618)

st.mcp <- hr_mcp(st.track, levels = c(0.5, 0.95))
st.mcp
head(st.mcp$data)
plot(st.mcp, col = c('red','green','blue'))
get_crs(st.mcp)

sf_sandtigerPositions <- sf::st_as_sf(sandtiger1, coords = c("longitude", "latitude"), 
                                    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% #"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "EPSG:4979"
  sf::st_transform(32618)


figure_sandtigerT1 <- #ggmap(SRWbasemap) +
  ggplot() +
  #minimum convex polygons
  geom_sf(data = st.mcp$mcp, 
          mapping = aes(fill = c(alpha("orange", 0.5), alpha("yellow", 0.2))), 
          size = 0.75, 
          linewidth = 0.9, 
          inherit.aes = FALSE) +
  #geom_sf(data = subset(sf_sandbarPositions, sf_sandbarPositions$dep_period), 
  #aes(color = "red"), alpha = 0.4, shape = 16) +
  #sandbar positions colored by year
  geom_sf(data = sf_sandtigerPositions, 
          mapping = aes(color = year), 
          alpha = 0.6, 
          shape = 16, 
          inherit.aes = FALSE) +
  #receiver stations showing target locations
  geom_sf(data = subset(sf_stations, sf_receivers$dep_period == "D1 target"), 
          mapping = aes(color = "black"), 
          size = 3.5, 
          inherit.aes = FALSE) +
  coord_sf(crs = st_crs(32618)) + #3857
  annotation_scale(location = "bl", width_hint = 0.13) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.45, "in"), pad_y = unit(0.01, "in"),
                         style = north_arrow_fancy_orienteering) +
  #coord_sf(xlim = c( st_bbox(sf_sandbarPositions)[["xmin"]]-400,  st_bbox(sf_sandbarPositions)[["xmax"]]+400), 
  #ylim = c( st_bbox(sf_sandbarPositions)[["ymin"]]-400,  st_bbox(sf_sandbarPositions)[["ymax"]]+400))  +
  scale_fill_manual(values = c(alpha("orange", 0.5), alpha("yellow", 0.2)), labels = c("50%", "95%")) + #these effectively make the legend
  scale_color_manual(values = c("red", "blue", "black"), 
                     labels = c("2022 positions", "2023 positions", "Target receiver locations")) + #these effectively make the legend
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
ggsave(paste0(owd,"/","figure_sandtigerT1.png"))

# sandbar MCP -------------------------------------------------------------------

sandbar1 <- orstedPosMerged2 %>% 
  filter(common_name_e == "Sandbar") %>% 
  mutate(est = as.POSIXct(est), 
         time = as.POSIXct(time),
         name = as.factor(name), 
         year = format(as.POSIXct(est), format = "%Y")) %>% 
  arrange(est) %>% 
  arrange(name)

sb.track <- make_track(sandbar1, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

str(sb.track)
class(sb.track)
sb.track
summary(sb.track)
get_crs(sb.track)
sb.track <- transform_coords(sb.track, crs_to=32618)

sb.mcp <- hr_mcp(sb.track, levels = c(0.5, 0.95))
sb.mcp
head(sb.mcp$data)
plot(sb.mcp, col = c('red','green','blue'))
get_crs(sb.mcp)

sf_sandbarPositions <- sf::st_as_sf(sandbar1, coords = c("longitude", "latitude"), 
                                    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% #"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "EPSG:4979"
  sf::st_transform(32618)


figure_sandbarT1 <- #ggmap(SRWbasemap) +
  ggplot() +
  #minimum convex polygons
  geom_sf(data = sb.mcp$mcp, 
          mapping = aes(fill = c(alpha("orange", 0.5), alpha("yellow", 0.2))), 
          size = 0.75, 
          linewidth = 0.9, 
          inherit.aes = FALSE) +
  #geom_sf(data = subset(sf_sandbarPositions, sf_sandbarPositions$dep_period), 
  #aes(color = "red"), alpha = 0.4, shape = 16) +
  #sandbar positions colored by year
  geom_sf(data = sf_sandbarPositions, 
          mapping = aes(color = year), 
          alpha = 0.6, 
          shape = 16, 
          inherit.aes = FALSE) +
  #receiver stations showing target locations
  geom_sf(data = subset(sf_stations, sf_receivers$dep_period == "D1 target"), 
          mapping = aes(color = "black"), 
          size = 3.5, 
          inherit.aes = FALSE) +
  coord_sf(crs = st_crs(32618)) + #3857
  annotation_scale(location = "bl", width_hint = 0.13) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.45, "in"), pad_y = unit(0.01, "in"),
                         style = north_arrow_fancy_orienteering) +
  #coord_sf(xlim = c( st_bbox(sf_sandbarPositions)[["xmin"]]-400,  st_bbox(sf_sandbarPositions)[["xmax"]]+400), 
  #ylim = c( st_bbox(sf_sandbarPositions)[["ymin"]]-400,  st_bbox(sf_sandbarPositions)[["ymax"]]+400))  +
  scale_fill_manual(values = c(alpha("orange", 0.5), alpha("yellow", 0.2)), labels = c("50%", "95%")) + #these effectively make the legend
  scale_color_manual(values = c("red", "blue", "black"), 
                     labels = c("2022 positions", "2023 positions", "Target receiver locations")) + #these effectively make the legend
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

# smooth dogfish MCP -------------------------------------------------------------------

smoothdog1 <- orstedPosMerged2 %>% 
  filter(common_name_e == "Smooth Dogfish") %>% 
  mutate(est = as.POSIXct(est), 
         time = as.POSIXct(time),
         name = as.factor(name), 
         year = format(as.POSIXct(est), format = "%Y")) %>% 
  arrange(est) %>% 
  arrange(name)

sd.track <- make_track(smoothdog1, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

str(sd.track)
class(sd.track)
sd.track
summary(sd.track)
get_crs(sd.track)
sd.track <- transform_coords(sd.track, crs_to=32618)

sd.mcp <- hr_mcp(sd.track, levels = c(0.5, 0.95))
sd.mcp
head(sd.mcp$data)
plot(sd.mcp, col = c('red','green','blue'))
get_crs(sd.mcp)

sf_smoothdogPositions <- sf::st_as_sf(smoothdog1, coords = c("longitude", "latitude"), 
                                  crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% #"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "EPSG:4979"
  sf::st_transform(32618)


figure_smoothdogT1 <- #ggmap(SRWbasemap) +
  ggplot() +
  #minimum convex polygons
  geom_sf(data = sd.mcp$mcp, 
          mapping = aes(fill = c(alpha("orange", 0.5), alpha("yellow", 0.2))), 
          size = 0.75, 
          linewidth = 0.9, 
          inherit.aes = FALSE) +
  #geom_sf(data = subset(sf_sandbarPositions, sf_sandbarPositions$dep_period), 
  #aes(color = "red"), alpha = 0.4, shape = 16) +
  #sandbar positions colored by year
  geom_sf(data = sf_smoothdogPositions, 
          mapping = aes(color = year), 
          alpha = 0.6, 
          shape = 16, 
          inherit.aes = FALSE) +
  #receiver stations showing target locations
  geom_sf(data = subset(sf_stations, sf_receivers$dep_period == "D1 target"), 
          mapping = aes(color = "black"), 
          size = 3.5, 
          inherit.aes = FALSE) +
  coord_sf(crs = st_crs(32618)) + #3857
  annotation_scale(location = "bl", width_hint = 0.13) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.45, "in"), pad_y = unit(0.01, "in"),
                         style = north_arrow_fancy_orienteering) +
  #coord_sf(xlim = c( st_bbox(sf_sandbarPositions)[["xmin"]]-400,  st_bbox(sf_sandbarPositions)[["xmax"]]+400), 
  #ylim = c( st_bbox(sf_sandbarPositions)[["ymin"]]-400,  st_bbox(sf_sandbarPositions)[["ymax"]]+400))  +
  scale_fill_manual(values = c(alpha("orange", 0.5), alpha("yellow", 0.2)), labels = c("50%", "95%")) + #these effectively make the legend
  scale_color_manual(values = c("red", "blue", "black"), 
                     labels = c("2022 positions", "2023 positions", "Target receiver locations")) + #these effectively make the legend
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white",
                                        color = "white"),
        plot.background = element_rect(fill = "white",
                                       color = "white")) +
  labs(title = expression(paste("Space use by smooth dogfish ", italic("Mustelus canis"))),
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
figure_smoothdogT1
ggsave(paste0(owd,"/","smoothdogfish_T1.png"))
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


# D1 missing receivers -------------------------------------------

############ working on this in untitled 3 (unless it gets accidentally deleted)
summary(as.factor(orstedPosMerged2$type))
syncTags1 <- orstedPosMerged2 %>% 
  filter(type == "Sync")
summary(as.factor(syncTags1$name))
summary(syncTags1$hp_em)

missingSyncTags_D1 <- orstedPosMerged2 %>% 
  filter(type == "Sync") %>% 
  filter(dep_period == "Deployment 1") %>% 
  filter(name == "R1C1" | name == "R1C2" | name == "R1C3" | name == "R1C4" | name == "R1C5" | name == "R1C6" | name == "R1C8" | name == "R2C6" |
           name == "R4C5" | name == "R4C7") %>% 
  mutate(time = as.POSIXct(time), 
         est = as.POSIXct(est))

firstLastD1 <- orstedPosMerged2 %>% 
  filter(type == "Sync") %>% 
  arrange(est) %>% 
  filter(dep_period == "Deployment 1") %>% 
  filter(name == "R1C1" | name == "R1C2" | name == "R1C3" | name == "R1C4" | name == "R1C5" | name == "R1C6" | name == "R1C8" | name == "R2C6" |
           name == "R4C5" | name == "R4C7") %>% 
  group_by(name) %>% 
  filter(est == min(est) | est == max(est)) %>% 
  arrange(est) %>% #gets everything in chronological order
  arrange(name) #gets the names in alphabetical order
firstLastD1$position <- rep(c("first position", "last position"), times = length(unique(firstLastD1$name)))

lostStationsPositionsD1 <- ggplot() +
  geom_point(data = missingSyncTags_D1, aes(x = longitude, y = latitude, color = name), alpha = 0.1) +
  geom_point(data = filter(stations, dep_period == "Deployment 1"), aes(x = longitude, y = latitude, shape = status), size = 3) +
  geom_point(data = firstLastD1, aes(x = longitude, y = latitude, fill = name, shape = position), size = 5) +
  scale_shape_manual(values = c("downloaded" = 19, "lost" = 4, "first position" = 21, "last position" = 22)) + 
  theme_minimal() +
  guides(colour = "none", #guide_legend(override.aes = list(alpha = 1)) 
         fill = "none") +
  theme(panel.background = element_rect(fill = "white",
                                        color = "white"),
        plot.background = element_rect(fill = "white",
                                       color = "white")) +
  labs(title = expression(paste("B)")),
       x =NULL, y = NULL, color = "stations", shape = "") + #, 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0), 
        axis.title.x = element_text(size= 14), 
        axis.text.x=element_text(size=12, color="black", hjust = 0.5),
        axis.title.y = element_text(size= 14),
        axis.text.y=element_text(size=12, color="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) +
  theme(legend.position = c(.89, .11),
        #legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
lostStationsPositionsD1

lostStationsAbacusD1 <- ggplot() +
  geom_point(data = missingSyncTags_D1, aes(x = time, y = name, color = name)) +
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
lostStationsAbacusD1

missingReceiversD1 <- ggarrange(lostStationsAbacusD1, lostStationsPositionsD1,
                              ncol = 2, nrow = 1)
missingReceiversD1
#annotate_figure(allSpeciesDetections, left = text_grob("Total detections", rot = 90, size = 20, face = "bold"),
#bottom = text_grob("Month", size = 20, face = "bold"))
ggsave(paste0(owd,"/","missingReceiversD1Positions.png"), width = 19, height = 10) 

# notes -------------------------------------------------------------------

#COA in VTrack
#igraph package for network analyses
#approach links the movements of each tagged individual to each visited area (Dakity or Manglar Bay) by an edge (arrow) and is weighted by the 
#number of movements/detections, thus providing an indicator of the use of each region
#produce bipartite graphs
#Brownscombe et al., () paper did this with bonefish comparing movements between bays, could use this for looking at movement between ACOE borrows
#state-space models using aniMotum R package (Jonsen et al., 2023)


