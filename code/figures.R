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
  labs(title = NULL, x =NULL, y = "HPEs")

plotRMSE <- orstedPosMerged2 %>% 
  filter(rmse < 10) %>% 
  group_by(dep_period) %>%
  ggplot() +
  geom_violin(aes(x = dep_period, y = rmse)) +
  theme_minimal() +
  labs(title = NULL, x =NULL, y = "RMSE")

errorComp <- ggarrange(plotHPEm, plotHPEs, plotRMSE, 
                       ncol = 3, nrow = 1)
errorComp

# figures ----------------------------------------------------------


# organizing station info -------------------------------------------------

stations <- read.csv("sunrise_coords.csv", header = TRUE) %>% 
  mutate(dep_period = as.factor(dep_period), 
         name = as.factor(name),
         status = as.factor(status))

sf_receivers <- sf::st_as_sf(stations, coords = c("longitude", "latitude"), 
                             crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% 
  #st_set_crs(32736) %>% 
  sf::st_transform(32618)

# deployment time periods -------------------------------------------------

D1start <- as.POSIXct("2022-07-30 12:00:00")
D1end <- as.POSIXct("2023-05-17 08:00:00")
D2start <- as.POSIXct("2023-05-17 15:55:00")
D2end <- as.POSIXct("2023-10-27 08:00:00")

# dusky MCP -------------------------------------------------------------------

#----- D1 -----#
duskyD1 <- orstedPosMerged2 %>% 
  filter(common_name_e == "Dusky") %>% 
  filter(est >= D1start & est < D1end)
str(duskyD1)

#test.track <- make_track(dusky1, longitude, latitude, time, crs=32618)
test.track <- make_track(duskyD1, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

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

sf_duskyPositions <- sf::st_as_sf(duskyD1, coords = c("longitude", "latitude"), 
                                  crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% #"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "EPSG:4979"
  sf::st_transform(32618)

#sf_duskyPositions <- sf::st_as_sf(dusky1, coords = c("longitude", "latitude")) %>% 
#sf::st_set_crs(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#sf::st_transform(32618)

confColor <- c(alpha("orange", 0.5), alpha("yellow", 0.2))
library(shades)
figure_duskyD1 <- ggplot() +
  geom_sf(data = dat.mcp$mcp, aes(fill = c(alpha("orange", 0.5), alpha("yellow", 0.2))), size = 0.75, linewidth = 0.9) +
  geom_sf(data = sf_duskyPositions, aes(color = "red"), alpha = 0.4, shape = 19) +
  geom_sf(data = subset(sf_receivers, sf_receivers$dep_period == "D1 target"), aes(color = "dark gray"), size = 2, shape = 17) +
  geom_sf(data = subset(sf_receivers, sf_receivers$dep_period == "Deployment 1"), aes(color = "black"), size = 3.5, shape = 19) + #this should be geom_sf but having CRS issues for some reason
  annotation_scale(location = "bl", width_hint = 0.13) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.45, "in"), pad_y = unit(0.01, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c( st_bbox(sf_duskyPositions)[["xmin"]]-0.001,  st_bbox(sf_duskyPositions)[["xmax"]]+0.001), 
           ylim = c( st_bbox(sf_duskyPositions)[["ymin"]]-0.001,  st_bbox(sf_duskyPositions)[["ymax"]]+0.001))  +
  theme_minimal() +
  scale_fill_manual(values = c(alpha("orange", 0.5), alpha("yellow", 0.2)), labels = c("50%", "95%")) + #these effectively make the legend
  #scale_color_manual(values = c("black", "dark gray", "red"), labels = c("Active Receivers", "Deployed Receivers", "Animal positions")) + 
  scale_colour_manual(name = "Receiver Status and Positions",
                      labels = c("Active Receivers", "Deployed Receivers", "Animal positions"),
                      values = c("black", "dark gray", "red")) +   
  scale_shape_manual(name = "Receiver Status and Positions",
                     labels = c("Active Receivers", "Deployed Receivers", "Animal positions"),
                     values = c(19, 17, 19)) +
  #these effectively make the legend
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white",
                                        color = "white"),
        plot.background = element_rect(fill = "white",
                                       color = "white")) +
  labs(title = "Deployment 1", x =NULL, y = NULL, fill = "", color = "") + #, 
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
figure_duskyD1
ggsave(paste0(owd,"/","duskysharks_D1.png"))

#----- D2 -----#
duskyD2 <- orstedPosMerged2 %>% 
  filter(common_name_e == "Dusky") %>% 
  filter(est >= D2start & est < D2end)
str(duskyD2)

#test.track <- make_track(dusky1, longitude, latitude, time, crs=32618)
test.track <- make_track(duskyD2, longitude, latitude, time, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

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

sf_duskyPositions <- sf::st_as_sf(duskyD2, coords = c("longitude", "latitude"), 
                                  crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% #"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "EPSG:4979"
  sf::st_transform(32618)

#sf_duskyPositions <- sf::st_as_sf(dusky1, coords = c("longitude", "latitude")) %>% 
#sf::st_set_crs(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#sf::st_transform(32618)

confColor <- c(alpha("orange", 0.5), alpha("yellow", 0.2))
library(shades)
figure_duskyD2 <- ggplot() +
  geom_sf(data = dat.mcp$mcp, aes(fill = c(alpha("orange", 0.5), alpha("yellow", 0.2))), size = 0.75, linewidth = 0.9) +
  geom_sf(data = sf_duskyPositions, aes(color = "red"), alpha = 0.4, shape = 19) +
  geom_sf(data = subset(sf_receivers, sf_receivers$dep_period == "D2 target"), aes(color = "dark gray"), size = 2, shape = 17) +
  geom_sf(data = subset(sf_receivers, sf_receivers$dep_period == "Deployment 2"), aes(color = "black"), size = 3.5, shape = 19) + #this should be geom_sf but having CRS issues for some reason
  annotation_scale(location = "bl", width_hint = 0.13) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.45, "in"), pad_y = unit(0.01, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c( st_bbox(sf_duskyPositions)[["xmin"]]-0.001,  st_bbox(sf_duskyPositions)[["xmax"]]+0.001), 
           ylim = c( st_bbox(sf_duskyPositions)[["ymin"]]-0.001,  st_bbox(sf_duskyPositions)[["ymax"]]+0.001))  +
  theme_minimal() +
  scale_fill_manual(values = c(alpha("orange", 0.5), alpha("yellow", 0.2)), labels = c("50%", "95%")) + #these effectively make the legend
  #scale_color_manual(values = c("black", "dark gray", "red"), labels = c("Active Receivers", "Deployed Receivers", "Animal positions")) + 
  scale_colour_manual(name = "Receiver Status and Positions",
                      labels = c("Active Receivers", "Deployed Receivers", "Animal positions"),
                      values = c("black", "dark gray", "red")) +   
  scale_shape_manual(name = "Receiver Status and Positions",
                     labels = c("Active Receivers", "Deployed Receivers", "Animal positions"),
                     values = c(19, 17, 19)) +
  #these effectively make the legend
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white",
                                        color = "white"),
        plot.background = element_rect(fill = "white",
                                       color = "white")) +
  labs(title = "Deployment 2", x =NULL, y = NULL, fill = "", color = "") + #, 
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
figure_duskyD2
ggsave(paste0(owd,"/","duskysharks_D2.png"))

duskyMCP <- ggarrange(figure_duskyD1, figure_duskyD2)
duskyMCP
annotate_figure(duskyMCP, top = text_grob(expression(paste("Space use by dusky sharks ", italic("Carcharhinus obscurus")))))
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
  arrange(est) %>% 
  arrange(name) 
firstLastD1$position <- rep(c("first", "last"), times = length(unique(firstLastD1$name)))

lostStationsPositionsD1 <- ggplot() +
  geom_point(data = missingSyncTags_D1, aes(x = longitude, y = latitude, color = name), alpha = 0.1) +
  geom_point(data = filter(stations, dep_period == "Deployment 1"), aes(x = longitude, y = latitude, shape = status), size = 3) +
  geom_point(data = firstLastD1, aes(x = longitude, y = latitude, fill = name, shape = position), size = 5) +
  scale_shape_manual(values = c("downloaded" = 19, "lost" = 4, "first" = 21, "last" = 22)) + 
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


# omitted -----------------------------------------------------------------

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

# deployment target station info
stations <- read_csv("sunrise_coords.csv") %>% 
  as.data.frame() %>% 
  dplyr::select(name, latitude, longitude, depth, dep_period) %>% 
  mutate(name = as.factor(name), 
         latitude = as.numeric(latitude),
         longitude = as.numeric(longitude), 
         depth = as.numeric(depth), 
         dep_period = as.factor(dep_period))
str(stations)

sf_receivers <- sf::st_as_sf(stations, coords = c("longitude", "latitude"), 
                             crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% 
  #st_set_crs(32736) %>% 
  sf::st_transform(32618)
mapview::mapview(sf_receivers)

test <- ggplot() +
  geom_sf(data = sf_receivers, 
          aes(color = "black"), size = 3.5) #this should be geom_sf but having CRS issues for some reason


lostStationsD1 <- c("R1C1", "R1C2", "R1C3", "R1C4", "R1C5", 
                    "R1C6", "R1C8", "R2C6", "R4C5", "R4C7")  
lostStationsD1.Lat <- c(40.731390, 40.732400, 40.733390, 40.734490, 40.735500, 
                        40.736500, 40.738590, 40.733510, 40.726530, 40.728580)
lostStationsD1.Long <- c(-72.842100, -72.839070, -72.836200, -72.833170, -72.830180, 
                         -72.827290, -72.821270, -72.826290, -72.824220, -72.821290)

lostStationCoords_D1 <- cbind(lostStationsD1, lostStationsD1.Lat, lostStationsD1.Long) %>% 
  as.data.frame() %>% 
  mutate(lostStationsD1.Lat = as.numeric(lostStationsD1.Lat), 
         lostStationsD1.Long = as.numeric(lostStationsD1.Long))

lostStationsPositionsD1 <- ggplot() +
  geom_point(data = missingSyncTags_D1, aes(x = longitude, y = latitude, color = name)) +
  #geom_point(data = stations, aes(x = Longitude, y = Latitude), color = "black") +
  geom_point(data = lostStationCoords_D1, aes(x = lostStationsD1.Long, y = lostStationsD1.Lat), color = "black", shape = 4, size = 4) +
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
lostStationsPositionsD1






