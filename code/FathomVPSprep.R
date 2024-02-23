#Orsted VPS data prep
#Author: Brianna Cahill
#Purpose: Cut down detection files to reflect receiver movement
#Input: detection_ORSTED2022_20240106.csv, detection_ORSTED2023_20240106.csv
#Output: more CSVs but cut down into bite sized pieces for VPS analyses

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
library(VTrack)
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

detections2022 <- read.csv ("detection_ORSTED2022_20240106.csv", header = TRUE)
detections2023 <- read.csv ("detection_ORSTED2023_20240106.csv", header = TRUE) #only until early november 2023
detections <- rbind(detections2022, detections2023)

summaryDets <- detections2022 %>%
  group_by(Station.Name, Receiver) %>%
  summarise(count = n())
  
  