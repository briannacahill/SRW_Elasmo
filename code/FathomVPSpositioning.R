#Orsted VPS data analysis
#Author: Brianna Cahill
#Purpose: Import animal position data within the SRW array and create figures
#Input: so many CSVs
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

# TO DO  -------------------------------------------------------------

  # gotta get the symbology figured out and how to show different shapes in the same legend
  # fix up sand tiger, sandbar and add smooth dogfish plots
  # lost receiver sync tag positions for D1 and D2, gotta fix this

# importing data ----------------------------------------------------------

  # this needs to be rerun EVERY time the positioning files are updated

#----- SRW deployment 1 -----#
orstedD1T1 <- list.files(path = "positioning/Orsted_D1T1_20221223",
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  lapply(\(x) mutate(x, across(SensorUnit, as.numeric))) %>% #some columns were logical while others were numeric, this allows them to be merged
  bind_rows() %>% 
  mutate(depPeriod = "Deployment 1") %>% 
  as.data.frame()
str(orstedD1T1)

orstedD1T2 <- list.files(path = "positioning/Orsted_D1T2_20230304",
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  lapply(\(x) mutate(x, across(SensorUnit, as.numeric))) %>% #some columns were logical while others were numeric, this allows them to be merged
  bind_rows() %>% 
  mutate(depPeriod = "Deployment 1") %>% 
  as.data.frame()
str(orstedD1T2)

#----- SRW deployment 2 -----#
orstedD2T3 <- list.files(path = "positioning/Orsted_D2T3_20230913",
                         pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  lapply(\(x) mutate(x, across(SensorUnit, as.numeric))) %>% #some columns were logical while others were numeric, this allows them to be merged
  bind_rows() %>% 
  mutate(depPeriod = "Deployment 2") %>% 
  as.data.frame()
str(orstedD2T3)

orstedAllPos <- rbind(orstedD1T1, orstedD1T2, 
                      orstedD2T3)

orstedAllPos$EST <- with_tz(orstedAllPos$Time, "America/New_York")

#filter out high error estimates HPEm (HPEm > 10; Bohaboy et al., 2022)
#orstedAllPos <- orstedAllPos %>% 
  #filter(HPEm < 10 | is.na(HPEm))
#summary(orstedAllPos$HPEm)

#need to skip first two empty rows, create column headers, then actually read in the CSV
  #doing this because I'm lazy and I don't feel like creating a CSV from something I already have to use in fathom position
  #probably easiest to use the most up to spec sheet containing tag information
myDeviceCols <- as.character(read_excel("FPspec_D1T1_20221223.xlsx", sheet = "Devices", skip = 2, n_max = 1, col_names = FALSE)) #need to skip first two empty
devices <- read_excel("FPspec_D1T1_20221223.xlsx", sheet = "Devices", skip = 2, col_names = myDeviceCols)
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

filteredPos <- orstedPosMerged2 %>% 
  filter(hp_em < 10) %>% #orstedPosMerged2
  as.data.frame()
str(filteredPos)

# summary and csv things --------------------------------------------------

#stony brook tag info
sbuFathomPositions <- orstedPosMerged2 %>% 
  filter(common_name_e == "Dusky" | common_name_e == "Sand Tiger" | common_name_e == "Sandbar" | common_name_e == "Smooth Dogfish" | common_name_e == "Winter Skate")
sbuFathomPositions
write.csv(sbuFathomPositions, paste0("sbuFathomPositions.csv")) 
write.csv(sbuFathomPositions, paste0(owd, "/", "sbuFathomPositions.csv")) 

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
write.csv(syncTagPositions, paste0("syncTagPositions.csv"))  
write.csv(syncTagPositions, paste0(owd, "/", "syncTagPositions.csv"))  

msRequest <- orstedPosMerged2 %>% 
  filter(full_id == "A69-1601-60134" | full_id == "A69-1601-60141" | full_id == "A69-1601-60219")
write.csv(msRequest, paste0(owd, "/", "msRequest_20240318.csv")) 
