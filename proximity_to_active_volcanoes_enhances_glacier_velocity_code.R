#### Script Info ####

# Statistical analyses for Mallalieu et al. 2023 "Proximity to active volcanoes enhances glacier velocity"
# Author: Dr Joseph Mallalieu, Research Fellow, University of Birmingham, UK
# Date: November 2023
# Written in: RStudio 2022.07.2+576 "Spotted Wakerobin" Release (e7373ef832b49b2a9b88162cfe7eac5f22c40b34, 2022-09-06) for Windows
# Available under a CC-BY-NC-SA license https://creativecommons.org/licenses/
# Contact: joe_mallalieu@outlook.com
# The data for these analyses can be downloaded from the UK Polar Data Centre: XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Notes on terminology: In the published manuscript the terms 'volcanic glacier' and 'non-volcanic glacier'...
# ...used here, are substituted for the terms 'glacier near a volcano' and 'other glacier' respectively.





#### Setup ####

# Load libraries for analysis
library("dplyr")
library("ggplot2")
library("lsr")
library("tidyr")
library("moments")
library("car")
library("caret")
library("corrplot")
library("ggpubr")
library("report")
library("ggpointdensity")
library("viridis")
library("scales")
library("lme4")
library("MuMIn")
library("MASS")
library("robustlmm")
library("nlme")

# Set scientific notation
options(scipen=999)   # removes scientific notation
options(scipen=0)   # reverts back to scientific notation

# Import data
volcglacdata <- read.csv(file = "put your file path here")
View(volcglacdata)





#### Data Management and Cleaning ####

# Add an absolute latitude variable (i.e. convert any negative latitude values to positive values)
volcglacdata$glacier.lat.abs <- abs(volcglacdata$glacier.lat)

# Change baseline of Eruption (aBP) variable to 2019 from 2018 (permits plotting of variable on log x axis, and removes problem of zero inflation in inferential stats)
volcglacdata$eruption.abp <- volcglacdata$eruption.abp + 1

# Convert variables to correct data type
volcglacdata$glacier.type <- as.factor(volcglacdata$glacier.type)
volcglacdata$glacier.no.volc.centres <- as.numeric(as.character(volcglacdata$glacier.no.volc.centres))
volcglacdata$glacier.distance.to.vent <- as.factor(volcglacdata$glacier.distance.to.vent)
volcglacdata$o1.region <- as.factor(volcglacdata$o1.region)
volcglacdata$o2.region <- as.factor(volcglacdata$o2.region)
volcglacdata$glacier.zmin <- as.numeric(as.character(volcglacdata$glacier.zmin))
volcglacdata$glacier.zmax <- as.numeric(as.character(volcglacdata$glacier.zmax))
volcglacdata$glacier.zmed <- as.numeric(as.character(volcglacdata$glacier.zmed))
volcglacdata$glacier.zrange <- as.numeric(as.character(volcglacdata$glacier.zrange))
volcglacdata$glacier.aspect <- as.numeric(as.character(volcglacdata$glacier.aspect))
volcglacdata$glacier.max.flowline <- as.numeric(as.character(volcglacdata$glacier.max.flowline))
volcglacdata$glacier.status <- as.factor(volcglacdata$glacier.status)
volcglacdata$glacier.connect <- as.factor(volcglacdata$glacier.connect)
volcglacdata$glacier.form <- as.factor(volcglacdata$glacier.form)
volcglacdata$glacier.surging <- as.factor(volcglacdata$glacier.surging)
volcglacdata$glacier.linkages <- as.factor(volcglacdata$glacier.linkages)
volcglacdata$glacier.thick.variety <- as.numeric(as.character(volcglacdata$glacier.thick.variety))
volcglacdata$volcano.number <- as.character(volcglacdata$volcano.number)
volcglacdata$volcano.country <- as.factor(volcglacdata$volcano.country)
volcglacdata$volcano.type <- as.factor(volcglacdata$volcano.type)
volcglacdata$eruption.evidence <- as.factor(volcglacdata$eruption.evidence)
volcglacdata$volcano.region <- as.factor(volcglacdata$volcano.region)
volcglacdata$volcano.subregion <- as.factor(volcglacdata$volcano.subregion)
volcglacdata$volcano.zmax <- as.numeric(as.character(volcglacdata$volcano.zmax))
volcglacdata$dominant.rock.type <- as.factor(volcglacdata$dominant.rock.type)
volcglacdata$tectonic.setting <- as.factor(volcglacdata$tectonic.setting)
volcglacdata$eruption.abp <- as.numeric(as.character(volcglacdata$eruption.abp))
volcglacdata$no.holocene.eruptions <- as.numeric(as.character(volcglacdata$no.holocene.eruptions))
volcglacdata$precip.1991.2020 <- as.numeric(as.character(volcglacdata$precip.1991.2020))
volcglacdata$precip.2017.2018 <- as.numeric(as.character(volcglacdata$precip.2017.2018))

# Consolidate number of factors in 'volcano.type' and eruption.evidence' variables
volcglacdata$volcano.type <- dplyr::recode(volcglacdata$volcano.type, "Pyroclastic cone(s)" = "Pyroclastic cone")   # merges the former into the latter
volcglacdata$volcano.type <- dplyr::recode(volcglacdata$volcano.type, "Shield(s)" = "Shield")
volcglacdata$volcano.type <- dplyr::recode(volcglacdata$volcano.type, "Stratovolcano(es)" = "Stratovolcano")
volcglacdata$volcano.type <- dplyr::recode(volcglacdata$volcano.type, "Stratovolcano?" = "Stratovolcano")
volcglacdata$volcano.type <- dplyr::recode(volcglacdata$volcano.type, "Compound" = "Complex")
volcglacdata$eruption.evidence <- dplyr::recode(volcglacdata$eruption.evidence, "Unrest / Holocene" = "Evidence Uncertain")

# Shorten names of 'dominant.rock.type' and 'tectonic.setting' factors
volcglacdata$dominant.rock.type <- dplyr::recode(volcglacdata$dominant.rock.type, "Andesite / Basaltic Andesite" = "Andesite/Bas.Andesite")
volcglacdata$dominant.rock.type <- dplyr::recode(volcglacdata$dominant.rock.type, "Basalt / Picro-Basalt" = "Basalt/Picro-Basalt")
volcglacdata$dominant.rock.type <- dplyr::recode(volcglacdata$dominant.rock.type, "Trachyandesite / Basaltic Trachyandesite" = "Trachyandesite/Bas.Trachyandesite")
volcglacdata$dominant.rock.type <- dplyr::recode(volcglacdata$dominant.rock.type, "Trachybasalt / Tephrite Basanite" = "Trachybasalt/Tephrite Basanite")
volcglacdata$dominant.rock.type <- dplyr::recode(volcglacdata$dominant.rock.type, "Trachyte / Trachydacite" = "Trachyte/Trachydacite")
volcglacdata$tectonic.setting <- dplyr::recode(volcglacdata$tectonic.setting, "Intraplate / Continental crust (>25 km)" = "Intraplate (cont. crust >25 km)")
volcglacdata$tectonic.setting <- dplyr::recode(volcglacdata$tectonic.setting, "Intraplate / Oceanic crust (< 15 km)" = "Intraplate (ocean. crust <15 km)")
volcglacdata$tectonic.setting <- dplyr::recode(volcglacdata$tectonic.setting, "Rift zone / Oceanic crust (< 15 km)" = "Rift Zone (ocean. crust <15 km)")
volcglacdata$tectonic.setting <- dplyr::recode(volcglacdata$tectonic.setting, "Subduction zone / Continental crust (>25 km)" = "Subduct. Zone (cont. crust >25 km)")
volcglacdata$tectonic.setting <- dplyr::recode(volcglacdata$tectonic.setting, "Subduction zone / Intermediate crust (15-25 km)" = "Subduct. Zone (intermed. crust 15-25 km)")
volcglacdata$tectonic.setting <- dplyr::recode(volcglacdata$tectonic.setting, "Subduction zone / Oceanic crust (< 15 km)" = "Subduct. Zone (ocean. crust <15 km)")

# Amend names of 'glacier.distance.to.vent' factors
volcglacdata$glacier.distance.to.vent <- dplyr::recode(volcglacdata$glacier.distance.to.vent, "> 1.0 km - <= 2.5 km" = "> 1.0 - <= 2.5 km")
volcglacdata$glacier.distance.to.vent <- dplyr::recode(volcglacdata$glacier.distance.to.vent, " > 2.5 km - <= 5km" = "> 2.5 - <= 5.0 km")

# Order factors by preference (determines reference factor in linear models)
volcglacdata$glacier.type <- factor(volcglacdata$glacier.type, levels = c("non-volcanic glacier", "volcanic glacier"))
volcglacdata$volcano.type <- factor(volcglacdata$volcano.type, levels = c("Stratovolcano", "Shield", "Complex", "Caldera", "Pyroclastic cone", "Subglacial", "Volcanic field"))
volcglacdata$eruption.evidence <- factor(volcglacdata$eruption.evidence, levels = c("Eruption Observed", "Eruption Dated", "Evidence Credible", "Evidence Uncertain"))
volcglacdata$tectonic.setting <- factor(volcglacdata$tectonic.setting, levels = c("Subduct. Zone (cont. crust >25 km)", "Subduct. Zone (intermed. crust 15-25 km)", "Rift Zone (ocean. crust <15 km)", "Intraplate (cont. crust >25 km)", "Intraplate (ocean. crust <15 km)"))
volcglacdata$glacier.distance.to.vent <- factor(volcglacdata$glacier.distance.to.vent, levels = c("<= 1.0 km", "> 1.0 - <= 2.5 km", "> 2.5 - <= 5.0 km"))

# Identify spurious glacier geometry values
volcglacdata %>% count(glacier.zmin < 0)
volcglacdata %>% count(glacier.zmax < 0)
volcglacdata %>% count(glacier.zmed < 0)
volcglacdata %>% count(glacier.zrange < 0)
volcglacdata %>% count(glacier.slope < 0)
volcglacdata %>% count(glacier.max.flowline < 0)
volcglacdata %>% count(glacier.aspect < 0)
volcglacdata %>% count(glacier.thick.median < 0.001)

# Set spurious glacier geometry values to NA
volcglacdata$glacier.zmin[which(volcglacdata$glacier.zmin < 0)] <- NA
volcglacdata$glacier.zmax[which(volcglacdata$glacier.zmax < 0)] <- NA
volcglacdata$glacier.zmed[which(volcglacdata$glacier.zmed < 0)] <- NA
volcglacdata$glacier.zrange[which(volcglacdata$glacier.zrange < 0)] <- NA
volcglacdata$glacier.slope[which(volcglacdata$glacier.slope < 0)] <- NA
volcglacdata$glacier.max.flowline[which(volcglacdata$glacier.max.flowline < 0)] <- NA
volcglacdata$glacier.aspect[which(volcglacdata$glacier.aspect < 0)] <- NA
volcglacdata$glacier.thick.median[which(volcglacdata$glacier.thick.median < 0.001)] <- NA

# Count number of glaciers with a velocity pixel count of <3
volcglacdata %>% filter(glacier.vel.count < 3) %>% summarise(n())
volcglacdata %>% group_by(glacier.type) %>% filter(glacier.vel.count < 3) %>% summarise(n())

# Count number of glaciers with a thickness pixel count of <3
volcglacdata %>% filter(glacier.thick.count < 3) %>% summarise(n())
volcglacdata %>% group_by(glacier.type) %>% filter(glacier.thick.count < 3) %>% summarise(n())

# Count number of glaciers with a velocity OR thickness pixel count of <3
volcglacdata %>% filter(glacier.vel.count < 3 | glacier.thick.count < 3) %>% summarise(n())
volcglacdata %>% group_by(glacier.type) %>% filter(glacier.vel.count < 3 | glacier.thick.count < 3) %>% summarise(n())

# Change all velocity and thickness data to NA where respective pixel count < 3
volcglacdata$glacier.vel.median[which(volcglacdata$glacier.vel.count < 3)] <- NA
volcglacdata$glacier.vel.sum[which(volcglacdata$glacier.vel.count < 3)] <- NA
volcglacdata$glacier.vel.mean[which(volcglacdata$glacier.vel.count < 3)] <- NA
volcglacdata$glacier.vel.stdev[which(volcglacdata$glacier.vel.count < 3)] <- NA
volcglacdata$glacier.vel.min[which(volcglacdata$glacier.vel.count < 3)] <- NA
volcglacdata$glacier.vel.max[which(volcglacdata$glacier.vel.count < 3)] <- NA
volcglacdata$glacier.vel.range[which(volcglacdata$glacier.vel.count < 3)] <- NA
volcglacdata$glacier.vel.minority[which(volcglacdata$glacier.vel.count < 3)] <- NA
volcglacdata$glacier.vel.majority[which(volcglacdata$glacier.vel.count < 3)] <- NA
volcglacdata$glacier.vel.variety[which(volcglacdata$glacier.vel.count < 3)] <- NA
volcglacdata$glacier.vel.variance[which(volcglacdata$glacier.vel.count < 3)] <- NA
volcglacdata$glacier.thick.median[which(volcglacdata$glacier.thick.count < 3)] <- NA
volcglacdata$glacier.thick.sum[which(volcglacdata$glacier.thick.count < 3)] <- NA
volcglacdata$glacier.thick.mean[which(volcglacdata$glacier.thick.count < 3)] <- NA
volcglacdata$glacier.thick.stdev[which(volcglacdata$glacier.thick.count < 3)] <- NA
volcglacdata$glacier.thick.min[which(volcglacdata$glacier.thick.count < 3)] <- NA
volcglacdata$glacier.thick.max[which(volcglacdata$glacier.thick.count < 3)] <- NA
volcglacdata$glacier.thick.range[which(volcglacdata$glacier.thick.count < 3)] <- NA
volcglacdata$glacier.thick.minority[which(volcglacdata$glacier.thick.count < 3)] <- NA
volcglacdata$glacier.thick.majority[which(volcglacdata$glacier.thick.count < 3)] <- NA
volcglacdata$glacier.thick.variety[which(volcglacdata$glacier.thick.count < 3)] <- NA
volcglacdata$glacier.thick.variance[which(volcglacdata$glacier.thick.count < 3)] <- NA

# Count number of glaciers with estimated geometry (ellipsoids) in the RGI (i.e. status > 0)
volcglacdata %>% filter(glacier.status != 0) %>% summarise(n())

# Change all glacier geometry, velocity, thickness and climate values to NA where glaciers have estimated geometry
volcglacdata$glacier.area[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.zmin[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.zmax[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.zmed[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.zrange[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.slope[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.aspect[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.max.flowline[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.vel.median[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.vel.sum[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.vel.mean[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.vel.stdev[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.vel.min[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.vel.max[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.vel.range[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.vel.minority[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.vel.majority[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.vel.variety[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.vel.variance[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.thick.median[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.thick.sum[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.thick.mean[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.thick.stdev[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.thick.min[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.thick.max[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.thick.range[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.thick.minority[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.thick.majority[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.thick.variety[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$glacier.thick.variance[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$temp.1991.2020[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$temp.2017.2018[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$precip.1991.2020[which(volcglacdata$glacier.status != 0)] <- NA
volcglacdata$precip.2017.2018[which(volcglacdata$glacier.status != 0)] <- NA

# Count number of glaciers strongly connected to the Greenland Ice Sheet (i.e. connection integer of > 1)
volcglacdata %>% filter(glacier.connect != 0, glacier.connect != 1) %>% summarise(n())

# Change all glacier geometry, velocity, thickness and climate values to NA where glaciers are strongly connected to the Greenland ice-sheet
volcglacdata$glacier.area[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.zmin[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.zmax[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.zmed[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.zrange[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.slope[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.aspect[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.max.flowline[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.vel.median[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.vel.sum[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.vel.mean[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.vel.stdev[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.vel.min[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.vel.max[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.vel.range[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.vel.minority[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.vel.majority[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.vel.variety[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.vel.variance[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.thick.median[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.thick.sum[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.thick.mean[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.thick.stdev[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.thick.min[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.thick.max[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.thick.range[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.thick.minority[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.thick.majority[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.thick.variety[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$glacier.thick.variance[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$temp.1991.2020[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$temp.2017.2018[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$precip.1991.2020[which(volcglacdata$glacier.connect == 2)] <- NA
volcglacdata$precip.2017.2018[which(volcglacdata$glacier.connect == 2)] <- NA

# Count number of glaciers with observed/probable/possible evidence of surging (i.e. surging integer of 1-3)
volcglacdata %>% filter(glacier.surging == 1 | glacier.surging == 2 | glacier.surging == 3 ) %>% summarise(n())

# Change all glacier geometry, velocity, thickness and climate values to NA where glaciers have demonstrated...
# ...observed/probable/possible evidence of surging
volcglacdata$glacier.area[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.zmin[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.zmax[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.zmed[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.zrange[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.slope[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.aspect[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.max.flowline[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.vel.median[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.vel.sum[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.vel.mean[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.vel.stdev[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.vel.min[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.vel.max[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.vel.range[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.vel.minority[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.vel.majority[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.vel.variety[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.vel.variance[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.thick.median[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.thick.sum[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.thick.mean[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.thick.stdev[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.thick.min[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.thick.max[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.thick.range[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.thick.minority[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.thick.majority[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.thick.variety[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$glacier.thick.variance[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$temp.1991.2020[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$temp.2017.2018[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$precip.1991.2020[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA
volcglacdata$precip.2017.2018[which(volcglacdata$glacier.surging == 1 | volcglacdata$glacier.surging == 2 | volcglacdata$glacier.surging == 3)] <- NA

# Count total number of glaciers in RGI dataset
length(volcglacdata$rgi.id)   ### n = 216876 # NB. 216876 - 489 (Ellipsoids) - 957 (connected to GrIS) - 1344 (surging) = 214086 in database

# Count total number of glaciers in RGI dataset by region
table.total.glac.region <- as.data.frame((volcglacdata %>% group_by(o1.region) %>% summarise(n())))
rownames(table.total.glac.region) <- c("Alaska", "W.Canada & USA", "Arctic Canada North", "Arctic Canada South",
                                       "Greenland Periphery", "Iceland", "Svalbard & Jan Mayen", "Scandinavia",
                                       "Russian Arctic", "North Asia", "Central Europe", "Caucasus & Middle East",
                                       "Central Asia", "South Asia West", "South Asia East", "Low Latitudes",
                                       "Southern Andes", "New Zealand", "Antarctic & Subantarctic")
print(table.total.glac.region)

# Count number of glacierised volcanoes in dataset (i.e. how many volcanoes have at least one glacier)
n_distinct(volcglacdata$volcano.number)

# Count number of non-volcanic and volcanic glaciers in RGI dataset
volcglacdata %>% group_by(glacier.type) %>% summarise(n())

# Table of non-volcanic and volcanic glaciers in RGI dataset by region
table.volcglac.count.region <- as.data.frame((volcglacdata %>% group_by(o1.region, glacier.type) %>% summarise(n())))
rownames(table.volcglac.count.region) <- c("Alaska", "Alaska (2)", "W.Canada & USA", "W.Canada & USA (2)", "Arctic Canada North",
                                           "Arctic Canada South", "Greenland Periphery", "Iceland", "Iceland (2)", "Svalbard & Jan Mayen",
                                           "Svalbard & Jan Mayen (2)", "Scandinavia", "Russian Arctic", "North Asia", "North Asia (2)",
                                           "Central Europe", "Caucasus & Middle East", "Caucasus & Middle East (2)", "Central Asia",
                                           "Central Asia (2)", "South Asia West", "South Asia East", "Low Latitudes", "Low Latitudes (2)",
                                           "Southern Andes", "Southern Andes (2)", "New Zealand", "New Zealand (2)",
                                           "Antarctic & Subantarctic", "Antarctic & Subantarctic (2)")
print(table.volcglac.count.region)

# Count total number of glaciers with median velocity values
sum(!is.na(volcglacdata$glacier.vel.median))

# Calculate % of glaciers that have median velocity value
((183546/216876)*100)

# Count total number of glaciers with median thickness value
sum(!is.na(volcglacdata$glacier.thick.median))

# Calculate % of glaciers that have median thickness value
((189736/216876)*100)





#### Comparison of 2017-2018 and 1991-2020 Climate ####

# Boxplot of mean annual air temperature in 1991-2020 and 2017-2018
boxplot(volcglacdata$temp.1991.2020, volcglacdata$temp.2017.2018,
        names=c("1991-2020", "2017-2018"),
        ylab = "Mean Annual Air Temperature (degrees)")

# Boxplot of total annual precipitation in 1991-2020 and 2017-2018
boxplot(volcglacdata$precip.1991.2020, volcglacdata$precip.2017.2018,
        names=c("1991-2020", "2017-2018"),
        ylab = "Total Annual Precipitation (mm)",
        log = "y")

# Calculate mean, median and standard deviation of temperature and precipitation for 1991-2020 and 2017-2018
t(volcglacdata %>% summarise(mean(na.omit(temp.1991.2020)),
                             mean(na.omit(temp.2017.2018)), median(na.omit(temp.1991.2020)), median(na.omit(temp.2017.2018)),
                             sd(na.omit(temp.1991.2020)), sd(na.omit(temp.2017.2018)), mean(na.omit(precip.1991.2020)),
                             mean(na.omit(precip.2017.2018)), median(na.omit(precip.1991.2020)), median(na.omit(precip.2017.2018)),
                             sd(na.omit(precip.1991.2020)), sd(na.omit(precip.2017.2018))))

# Calculate Cohen's D to investigate effect size between mean temperature in 1992-2020 and 2017-2018
cohensD(volcglacdata$temp.1991.2020, volcglacdata$temp.2017.2018)

# Calculate Cohen's D to investigate effect size between mean precipitation in 1992-2020 and 2017-2018
cohensD(volcglacdata$precip.1991.2020, volcglacdata$precip.2017.2018)





#### Non-Volcanic/Volcanic Glacier: Descriptive Statistics ####

# Summary statistics for all glaciers in dataset: MEANS
table.average.glacier <- t(as.data.frame(volcglacdata %>% summarise(n(),
                                                                    mean(glacier.long), mean(glacier.lat), mean(glacier.lat.abs), mean(glacier.area, na.rm = TRUE), mean(glacier.zmin, na.rm = TRUE),
                                                                    mean(glacier.zmax, na.rm = TRUE), mean(glacier.zmed, na.rm = TRUE), mean(glacier.zrange, na.rm = TRUE),
                                                                    mean(glacier.slope, na.rm = TRUE), mean(glacier.aspect, na.rm = TRUE), mean(glacier.max.flowline, na.rm = TRUE),
                                                                    mean(glacier.vel.count, na.rm = TRUE), mean(glacier.vel.mean, na.rm = TRUE),
                                                                    mean(glacier.vel.median, na.rm = TRUE), mean(glacier.thick.count, na.rm = TRUE), 
                                                                    mean(glacier.thick.mean, na.rm = TRUE), mean(glacier.thick.median, na.rm = TRUE), 
                                                                    mean(temp.2017.2018, na.rm = TRUE), mean(precip.2017.2018, na.rm = TRUE))
                                         %>% mutate_if(is.numeric, ~round(., 2))))
print(table.average.glacier)

# Summary statistics for all glaciers in dataset: MEDIANS
table.average.glacier.median <- t(as.data.frame(volcglacdata %>% summarise(n(),
                                                                           median(glacier.long), median(glacier.lat), median(glacier.lat.abs), median(glacier.area, na.rm = TRUE), median(glacier.zmin, na.rm = TRUE),
                                                                           median(glacier.zmax, na.rm = TRUE), median(glacier.zmed, na.rm = TRUE), median(glacier.zrange, na.rm = TRUE),
                                                                           median(glacier.slope, na.rm = TRUE), median(glacier.aspect, na.rm = TRUE), median(glacier.max.flowline, na.rm = TRUE),
                                                                           median(glacier.vel.count, na.rm = TRUE), median(glacier.vel.median, na.rm = TRUE),
                                                                           median(glacier.vel.median, na.rm = TRUE), median(glacier.thick.count, na.rm = TRUE), 
                                                                           median(glacier.thick.median, na.rm = TRUE), median(glacier.thick.median, na.rm = TRUE), 
                                                                           median(temp.2017.2018, na.rm = TRUE), median(precip.2017.2018, na.rm = TRUE))
                                                %>% mutate_if(is.numeric, ~round(., 2))))
print(table.average.glacier.median)

# Summary statistics of all glaciers by glacier type (non-volcanic/volcanic): MEANS
table.average.glacier.by.type <- t(as.data.frame(volcglacdata %>% group_by(glacier.type) %>% summarise(n(),
                                                                                                       mean(glacier.long), mean(glacier.lat), mean(glacier.lat.abs), mean(glacier.area, na.rm = TRUE), mean(glacier.zmin, na.rm = TRUE),
                                                                                                       mean(glacier.zmax, na.rm = TRUE), mean(glacier.zmed, na.rm = TRUE), mean(glacier.zrange, na.rm = TRUE),
                                                                                                       mean(glacier.slope, na.rm = TRUE), mean(glacier.aspect, na.rm = TRUE), mean(glacier.max.flowline, na.rm = TRUE),
                                                                                                       mean(glacier.vel.count, na.rm = TRUE), mean(glacier.vel.mean, na.rm = TRUE),
                                                                                                       mean(glacier.vel.median, na.rm = TRUE), mean(glacier.thick.count, na.rm = TRUE), mean(glacier.thick.mean, na.rm = TRUE),
                                                                                                       mean(glacier.thick.median, na.rm = TRUE), mean(temp.2017.2018, na.rm = TRUE), mean(precip.2017.2018, na.rm = TRUE)))
                                   %>% mutate_if(is.numeric, ~round(., 2)))
print(table.average.glacier.by.type)

# Summary statistics of all glaciers by glacier type (non-volcanic/volcanic): MEDIANS
table.average.glacier.by.type.median <- t(as.data.frame(volcglacdata %>% group_by(glacier.type) %>% summarise(n(),
                                                                                                              median(glacier.long), median(glacier.lat), median(glacier.lat.abs), median(glacier.area, na.rm = TRUE), median(glacier.zmin, na.rm = TRUE),
                                                                                                              median(glacier.zmax, na.rm = TRUE), median(glacier.zmed, na.rm = TRUE), median(glacier.zrange, na.rm = TRUE),
                                                                                                              median(glacier.slope, na.rm = TRUE), median(glacier.aspect, na.rm = TRUE), median(glacier.max.flowline, na.rm = TRUE),
                                                                                                              median(glacier.vel.count, na.rm = TRUE), median(glacier.vel.median, na.rm = TRUE),
                                                                                                              median(glacier.vel.median, na.rm = TRUE), median(glacier.thick.count, na.rm = TRUE),
                                                                                                              median(glacier.thick.median, na.rm = TRUE), median(glacier.thick.median, na.rm = TRUE),
                                                                                                              median(temp.2017.2018, na.rm = TRUE), median(precip.2017.2018, na.rm = TRUE)))
                                          %>% mutate_if(is.numeric, ~round(., 2)))
print(table.average.glacier.by.type.median)

# Calculate Cohen's D to investigate effect sizes between non-volcanic glacier and volcanic glacier variables
cohend.nonvolcglac <- volcglacdata %>% filter(glacier.type == "non-volcanic glacier")
cohend.volcglac <- volcglacdata %>% filter(glacier.type == "volcanic glacier")
cohensD(cohend.nonvolcglac$glacier.lat.abs, cohend.volcglac$glacier.lat.abs)
cohensD(cohend.nonvolcglac$glacier.area,cohend.volcglac$glacier.area)
cohensD(cohend.nonvolcglac$glacier.zmed,cohend.volcglac$glacier.zmed)
cohensD(cohend.nonvolcglac$glacier.slope,cohend.volcglac$glacier.slope)
cohensD(cohend.nonvolcglac$glacier.max.flowline,cohend.volcglac$glacier.max.flowline)
cohensD(cohend.nonvolcglac$glacier.vel.median, cohend.volcglac$glacier.vel.median)
cohensD(cohend.nonvolcglac$glacier.vel.stdev,cohend.volcglac$glacier.vel.stdev)
cohensD(cohend.nonvolcglac$glacier.thick.median,cohend.volcglac$glacier.thick.median)
cohensD(cohend.nonvolcglac$glacier.thick.stdev,cohend.volcglac$glacier.thick.stdev)
cohensD(cohend.nonvolcglac$temp.2017.2018,cohend.volcglac$temp.2017.2018)
cohensD(cohend.nonvolcglac$precip.2017.2018,cohend.volcglac$precip.2017.2018)

# Violin plots of non-volcanic and volcanic glaciers by variable (including means added as points). NB. 'Non-finite values' in warning messages are just NAs
# Create function for plotting minor log breaks
minor_breaks_log <- function(base) {
  force(base) 
  function(limits) {
    ggplot2:::calc_logticks(
      base = base, 
      minpow = floor(log(limits[1], base = base)), 
      maxpow = ceiling(log(limits[2], base = base))
    )$value
  }
}

# Latitude
mean.nonvolc.lat <- mean(volcglacdata$glacier.lat[volcglacdata$glacier.type == "non-volcanic glacier"], na.rm = TRUE)
mean.volc.lat <- mean(volcglacdata$glacier.lat[volcglacdata$glacier.type == "volcanic glacier"], na.rm = TRUE)
n.nonvolc.lat <- volcglacdata %>% filter(glacier.type == "non-volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.lat)))
n.volc.lat <- volcglacdata %>% filter(glacier.type == "volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.lat)))
violin.lat <- ggplot(volcglacdata, aes(x = glacier.type, y = glacier.lat, fill = glacier.type)) +
  geom_violin() +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "", y = "Latitude (°)") +
  scale_x_discrete(labels = c("Glaciers near a volcano", "Other glaciers")) +
  guides(fill = "none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Pastel1", direction = 1) +
  geom_point(aes(x = "non-volcanic glacier", y = mean.nonvolc.lat), col = "black", shape = 18, size = 4) +
  geom_point(aes(x = "volcanic glacier", y = mean.volc.lat), col = "black", shape = 18, size = 4)
violin.lat

# Absolute latitude
mean.nonvolc.abslat <- mean(volcglacdata$glacier.lat.abs[volcglacdata$glacier.type == "non-volcanic glacier"], na.rm = TRUE)
mean.volc.abslat <- mean(volcglacdata$glacier.lat.abs[volcglacdata$glacier.type == "volcanic glacier"], na.rm = TRUE)
n.nonvolc.abslat <- volcglacdata %>% filter(glacier.type == "non-volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.lat.abs)))
n.volc.abslat <- volcglacdata %>% filter(glacier.type == "volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.lat.abs)))
violin.abslat <- ggplot(volcglacdata, aes(x = glacier.type, y = glacier.lat.abs, fill = glacier.type)) +
  geom_violin() +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "", y = "Absolute Latitude (°)") +
  scale_x_discrete(labels = c("Glaciers near a volcano", "Other glaciers")) +
  guides(fill = "none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Pastel1", direction = 1) +
  geom_point(aes(x = "non-volcanic glacier", y = mean.nonvolc.abslat), col = "black", shape = 18, size = 4) +
  geom_point(aes(x = "volcanic glacier", y = mean.volc.abslat), col = "black", shape = 18, size = 4)
violin.abslat

# Area
mean.nonvolc.area <- mean(volcglacdata$glacier.area[volcglacdata$glacier.type == "non-volcanic glacier"], na.rm = TRUE)
mean.volc.area <- mean(volcglacdata$glacier.area[volcglacdata$glacier.type == "volcanic glacier"], na.rm = TRUE)
n.nonvolc.area <- volcglacdata %>% filter(glacier.type == "non-volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.area)))
n.volc.area <- volcglacdata %>% filter(glacier.type == "volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.area)))
violin.area <- ggplot(volcglacdata, aes(x=glacier.type, y=glacier.area, fill=glacier.type)) +
  geom_violin() +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "", y = expression(Area~~(km^2))) +
  scale_x_discrete(labels = c("Glaciers near a volcano", "Other glaciers")) +
  guides(fill="none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Pastel1", direction=1) +
  geom_point(aes(x="non-volcanic glacier",y=mean.nonvolc.area), col="black", shape = 18, size = 4) +
  geom_point(aes(x="volcanic glacier",y=mean.volc.area), col="black", shape = 18, size = 4)
violin.area

# Median elevation
mean.nonvolc.zmed <- mean(volcglacdata$glacier.zmed[volcglacdata$glacier.type == "non-volcanic glacier"], na.rm = TRUE)
mean.volc.zmed <- mean(volcglacdata$glacier.zmed[volcglacdata$glacier.type == "volcanic glacier"], na.rm = TRUE)
n.nonvolc.zmed <- volcglacdata %>% filter(glacier.type == "non-volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.zmed)))
n.volc.zmed <- volcglacdata %>% filter(glacier.type == "volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.zmed)))
violin.zmed <- ggplot(volcglacdata, aes(x = glacier.type, y = glacier.zmed, fill = glacier.type)) +
  geom_violin() +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "", y = "Median elevation (m)") +
  scale_x_discrete(labels = c("Glaciers near a volcano", "Other glaciers")) +
  guides(fill = "none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Pastel1", direction = 1) +
  geom_point(aes(x = "non-volcanic glacier", y = mean.nonvolc.zmed), col = "black", shape = 18, size = 4) +
  geom_point(aes(x = "volcanic glacier", y = mean.volc.zmed), col = "black", shape = 18, size = 4)
violin.zmed

# Slope
mean.nonvolc.slope <- mean(volcglacdata$glacier.slope[volcglacdata$glacier.type == "non-volcanic glacier"], na.rm = TRUE)
mean.volc.slope <- mean(volcglacdata$glacier.slope[volcglacdata$glacier.type == "volcanic glacier"], na.rm = TRUE)
n.nonvolc.slope <- volcglacdata %>% filter(glacier.type == "non-volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.slope)))
n.volc.slope <- volcglacdata %>% filter(glacier.type == "volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.slope)))
violin.slope <- ggplot(volcglacdata, aes(x = glacier.type, y = glacier.slope, fill = glacier.type)) +
  geom_violin() +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "", y = "Slope (°)") +
  scale_x_discrete(labels = c("Glaciers near a volcano", "Other glaciers")) +
  guides(fill = "none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Pastel1", direction = 1) +
  geom_point(aes(x = "non-volcanic glacier", y = mean.nonvolc.slope), col = "black", shape = 18, size = 4) +
  geom_point(aes(x = "volcanic glacier", y = mean.volc.slope), col = "black", shape = 18, size = 4)
violin.slope

# Maximum surface flowline length
mean.nonvolc.max.flowline <- mean(volcglacdata$glacier.max.flowline[volcglacdata$glacier.type == "non-volcanic glacier"], na.rm = TRUE)
mean.volc.max.flowline <- mean(volcglacdata$glacier.max.flowline[volcglacdata$glacier.type == "volcanic glacier"], na.rm = TRUE)
n.nonvolc.max.flowline <- volcglacdata %>% filter(glacier.type == "non-volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.max.flowline)))
n.volc.max.flowline <- volcglacdata %>% filter(glacier.type == "volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.max.flowline)))
violin.max.flowline <- ggplot(volcglacdata, aes(x=glacier.type, y=glacier.max.flowline, fill=glacier.type)) +
  geom_violin() +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "", y = "Maximum surface\n flowline length (m)") +
  scale_x_discrete(labels = c("Glaciers near a volcano", "Other glaciers")) +
  guides(fill="none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Pastel1", direction = 1) +
  geom_point(aes(x="non-volcanic glacier",y=mean.nonvolc.max.flowline), col="black", shape = 18, size = 4) +
  geom_point(aes(x="volcanic glacier",y=mean.volc.max.flowline), col="black", shape = 18, size = 4)
violin.max.flowline

# Median velocity
mean.nonvolc.vel.median <- mean(volcglacdata$glacier.vel.median[volcglacdata$glacier.type == "non-volcanic glacier"], na.rm = TRUE)
mean.volc.vel.median <- mean(volcglacdata$glacier.vel.median[volcglacdata$glacier.type == "volcanic glacier"], na.rm = TRUE)
n.nonvolc.vel.median <- volcglacdata %>% filter(glacier.type == "non-volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.vel.median)))
n.volc.vel.median <- volcglacdata %>% filter(glacier.type == "volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.vel.median)))
violin.vel.median <- ggplot(volcglacdata, aes(x=glacier.type, y=glacier.vel.median, fill=glacier.type)) +
  geom_violin() +
  scale_y_continuous(trans = "log10", labels = scales::comma, minor_breaks = minor_breaks_log(10)) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "", y = expression(Median~velocity~~(ma^-1))) +
  scale_x_discrete(labels = c("Glaciers near a volcano", "Other glaciers")) +
  guides(fill="none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Pastel1", direction = 1) +
  geom_point(aes(x="non-volcanic glacier",y=mean.nonvolc.vel.median), col="black", shape = 18, size = 4) +
  geom_point(aes(x="volcanic glacier",y=mean.volc.vel.median), col="black", shape = 18, size = 4)
violin.vel.median

### Median thickness. Note y axis cropped to 300 m...
volcglacdata %>% count(glacier.thick.median > 300)
mean.nonvolc.thick.median <- mean(volcglacdata$glacier.thick.median[volcglacdata$glacier.type == "non-volcanic glacier"], na.rm = TRUE)
mean.volc.thick.median <- mean(volcglacdata$glacier.thick.median[volcglacdata$glacier.type == "volcanic glacier"], na.rm = TRUE)
n.nonvolc.thick.median <- volcglacdata %>% filter(glacier.type == "non-volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.thick.median)))
n.volc.thick.median <- volcglacdata %>% filter(glacier.type == "volcanic glacier") %>% summarise(total_non_na = sum(!is.na(glacier.thick.median)))
violin.thick.median.nonlog <- ggplot(volcglacdata, aes(x = glacier.type, y = glacier.thick.median, fill = glacier.type)) +
  geom_violin() +
  ylim(c(0,300)) +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "", y = "Median thickness (m)") +
  scale_x_discrete(labels = c("Glaciers near a volcano", "Other glaciers")) +
  guides(fill = "none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Pastel1", direction = 1) +
  geom_point(aes(x = "non-volcanic glacier", y = mean.nonvolc.thick.median), col = "black", shape = 18, size = 4) +
  geom_point(aes(x = "volcanic glacier", y = mean.volc.thick.median), col = "black", shape = 18, size = 4)
violin.thick.median.nonlog

# Temperature
mean.nonvolc.temp <- mean(volcglacdata$temp.2017.2018[volcglacdata$glacier.type == "non-volcanic glacier"], na.rm = TRUE)
mean.volc.temp <- mean(volcglacdata$temp.2017.2018[volcglacdata$glacier.type == "volcanic glacier"], na.rm = TRUE)
n.nonvolc.temp <- volcglacdata %>% filter(glacier.type == "non-volcanic glacier") %>% summarise(total_non_na = sum(!is.na(temp.2017.2018)))
n.volc.temp <- volcglacdata %>% filter(glacier.type == "volcanic glacier") %>% summarise(total_non_na = sum(!is.na(temp.2017.2018)))
violin.temp <- ggplot(volcglacdata, aes(x = glacier.type, y = temp.2017.2018, fill = glacier.type)) +
  geom_violin() +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "", y = "Mean annual temperature\n 2017-2018 (°C)") +
  scale_x_discrete(labels = c("Glaciers near a volcano", "Other glaciers")) +
  guides(fill = "none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Pastel1", direction = 1) +
  geom_point(aes(x = "non-volcanic glacier", y = mean.nonvolc.temp), col = "black", shape = 18, size = 4) +
  geom_point(aes(x = "volcanic glacier", y = mean.volc.temp), col = "black", shape = 18, size = 4)
violin.temp

# Precipitation
mean.nonvolc.precip <- mean(volcglacdata$precip.2017.2018[volcglacdata$glacier.type == "non-volcanic glacier"], na.rm = TRUE)
mean.volc.precip <- mean(volcglacdata$precip.2017.2018[volcglacdata$glacier.type == "volcanic glacier"], na.rm = TRUE)
n.nonvolc.precip <- volcglacdata %>% filter(glacier.type == "non-volcanic glacier") %>% summarise(total_non_na = sum(!is.na(precip.2017.2018)))
n.volc.precip <- volcglacdata %>% filter(glacier.type == "volcanic glacier") %>% summarise(total_non_na = sum(!is.na(precip.2017.2018)))
violin.precip <- ggplot(volcglacdata, aes(x=glacier.type, y=precip.2017.2018, fill=glacier.type)) +
  geom_violin() +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "", y = "Total annual precipitation\n 2017-2018 (mm)") +
  scale_x_discrete(labels = c("Glaciers near a volcano", "Other glaciers")) +
  guides(fill="none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Pastel1", direction=1) +
  geom_point(aes(x="non-volcanic glacier",y=mean.nonvolc.precip), col="black", shape = 18, size = 4) +
  geom_point(aes(x="volcanic glacier",y=mean.volc.precip), col="black", shape = 18, size = 4)
violin.precip

# Composite plot of above violin plots
ggarrange(violin.abslat, violin.zmed, violin.area, violin.max.flowline, violin.thick.median.nonlog, 
          violin.slope, violin.temp, violin.precip,
          ncol = 2, nrow = 4, align = "hv")

# Box plot of median glacier velocity for all regions with non-volcanic and volcanic glaciers (excluding NZ due to absence of data). NB Outliers removed.
boxplot.region.vel <- volcglacdata %>% filter(o1.region == 1 | o1.region == 2 | o1.region == 6 | o1.region == 7 | o1.region == 10 | o1.region == 12 | o1.region == 13 | o1.region == 16 | o1.region == 17 | o1.region == 19)
boxplot.region.vel.plot <- ggplot(data = (boxplot.region.vel), aes(x=factor(o1.region), y=glacier.vel.median, fill=factor(glacier.type))) + 
  geom_boxplot(colour = "grey25", outlier.shape = NA) +
  scale_y_continuous(trans = "log10", labels = scales::comma, minor_breaks = minor_breaks_log(10), limits = c(0.5,200)) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  labs(x = "", y = expression(Median~velocity~~(ma^-1))) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_blank(),
        legend.position = ("top"), legend.title = element_blank(), legend.box.margin = margin(0, 0, -35, 0),
        legend.text = element_text(size = 11)) +
  scale_fill_brewer(palette = "Pastel1", direction=1, labels = c("Glaciers near a volcano", "Other glaciers")) +
  scale_x_discrete(guide = guide_axis(n.dodge=2),
                   labels=c("Alaska", "W.Can. & USA", "Iceland", "Svalbard & J.Mayen",
                            "North Asia", "Cauc. & Mid. East", "Central Asia", "Low Latitudes",
                            "Southern Andes", "Ant. & Subant."))
boxplot.region.vel.plot

# Composite plot of above velocity and regional velocity plots
ggarrange(violin.vel.median, boxplot.region.vel.plot, ncol = 1, align = "v")





#### Non-Volcanic/Volcanic Glacier: Inferential Statistics ####

# Create new data frame without glaciers that return NA for median glacier velocity
volcglacdata.inf <- volcglacdata %>% drop_na(glacier.vel.median)

# Histogram of dependent variable to determine distribution
hist(volcglacdata.inf$glacier.vel.median, breaks = 100)   # I.e. dependent variable possesses an inverse gaussian distribution
skewness(volcglacdata.inf$glacier.vel.median)
qqPlot(volcglacdata.inf$glacier.vel.median)

# Apply log10 transformation to correct for inverse gaussian distribution
hist(log10(volcglacdata.inf$glacier.vel.median), breaks = 100)
skewness(log10(volcglacdata.inf$glacier.vel.median))
qqPlot(log10(volcglacdata.inf$glacier.vel.median))

# Check for sparse variables in volcglacdata.inf dataframe
nearZeroVar(volcglacdata.inf)
colnames(volcglacdata.inf)

# Create new dataframe to check for multicollinearity in potential model variables
volcglacdata.inf.cor <- volcglacdata.inf

# Convert factor variables to numeric variables for correlation matrix
volcglacdata.inf.cor$glacier.type.integer <- as.numeric(as.character(volcglacdata.inf.cor$glacier.type.integer))
volcglacdata.inf.cor$o1.region <- as.numeric(as.character(volcglacdata.inf.cor$o1.region))
volcglacdata.inf.cor$glacier.status <- as.numeric(as.character(volcglacdata.inf.cor$glacier.status))
volcglacdata.inf.cor$glacier.connect <- as.numeric(as.character(volcglacdata.inf.cor$glacier.connect))
volcglacdata.inf.cor$glacier.form <- as.numeric(as.character(volcglacdata.inf.cor$glacier.form))
volcglacdata.inf.cor$glacier.surging <- as.numeric(as.character(volcglacdata.inf.cor$glacier.surging))
volcglacdata.inf.cor$glacier.linkages <- as.numeric(as.character(volcglacdata.inf.cor$glacier.linkages))

# Correlation matrix for all potential model variables:
correlation.potential.variables <- cor(volcglacdata.inf.cor[c("glacier.type.integer", "glacier.long", "glacier.lat",
                                                              "glacier.lat.abs", "o1.region", "glacier.area",
                                                              "glacier.zmin", "glacier.zmax", "glacier.zmed", "glacier.zrange", "glacier.slope",
                                                              "glacier.aspect", "glacier.max.flowline", "glacier.connect", "glacier.form", 
                                                              "glacier.surging", "glacier.linkages", "glacier.vel.median",
                                                              "glacier.thick.median", "temp.2017.2018", "precip.2017.2018")], use = "pairwise.complete.obs")
LM.corrplot.a <- corrplot(cor(correlation.potential.variables), method = "color", diag = FALSE, outline = TRUE, type = "lower", tl.col = "black", addCoef.col = TRUE, col = COL2('RdBu', 100), tl.srt = 90, cl.ratio = 0.2, tl.cex = 0.7, number.cex = 0.7)

# Correlation matrix for likely model variables:
correlation.likely.variables <- cor(volcglacdata.inf.cor[c("glacier.type.integer", "glacier.lat.abs", "glacier.zmed",
                                                           "glacier.slope", "glacier.vel.median",
                                                           "glacier.thick.median", "temp.2017.2018", "precip.2017.2018")], use = "pairwise.complete.obs")
LM.corrplot.b <- cor(correlation.likely.variables)
colnames(LM.corrplot.b) <- c("Glacier type", "Absolute latitude", "Median elevation", "Slope", "Median velocity", "Median thickness", "Mean annual\n temp. (2017-2018)", "Total annual\n precip. (2017-18)")
rownames(LM.corrplot.b) <- c("Glacier type", "Absolute latitude", "Median elevation", "Slope", "Median velocity", "Median thickness", "Mean annual\n temp. (2017-2018)", "Total annual\n precip. (2017-18)")
corrplot(LM.corrplot.b, method = "color", diag = FALSE, outline = TRUE, type = "lower", tl.col = "black", addCoef.col = TRUE, col = COL2('RdBu', 100), tl.srt = 90, cl.ratio = 0.2, tl.cex = 0.8, number.cex = 0.8)

# Investigate distribution of temperature and precipitation variables ahead of Principal Component Analysis
hist(volcglacdata.inf$temp.2017.2018)
skewness(volcglacdata.inf$temp.2017.2018)
hist(volcglacdata.inf$precip.2017.2018)
skewness(volcglacdata.inf$precip.2017.2018)

# Apply log transformation to correct skew in precipitation data
hist(log10(volcglacdata.inf$precip.2017.2018))
skewness(log10(volcglacdata.inf$precip.2017.2018))

# Create climate principal component (PC)
climatepc <- prcomp(~ temp.2017.2018 + log10(precip.2017.2018),
                    data = volcglacdata.inf,
                    scale = TRUE)
summary(climatepc)

# Calculate PC eigenvalues (only useful if eigenvalue > 1)
plot(climatepc, type = "lines")
climatepc$sdev^2

# Calculate PC eigenvectors (relative contribution of each variable to each principal component)
climatepc$rotation
biplot(climatepc, cex = 0.5) # visualises the loadings (the contribution of each variable to the PC)
biplot(climatepc, cex = 0.5, choices = c(1,2)) # plots PC1 versus PC2

# Define PC1 as a new variable
climate.pca.2017.2018 <- predict(climatepc)[,1]

# Add climate PC to dataframe
volcglacdata.inf <- cbind(volcglacdata.inf, climate.pca.2017.2018)

# Variance Inflation Factor (VIF) to check for multicollinearity between likely numerical and factor variables
lm.vif <- lm(glacier.vel.median ~ glacier.type + glacier.zmed + glacier.slope + 
               glacier.thick.median + climate.pca.2017.2018,  data = volcglacdata.inf)
as.data.frame(vif(lm.vif)) # All values <5, so no evidence of multicollinearity between factor variables and numerical variables

# Subset dataframe with likely model variables that require scaling and centring
volcglacdata.inf.sub <- subset(volcglacdata.inf, select = c("glacier.zmed", "glacier.slope", "glacier.thick.median",
                                                            "climate.pca.2017.2018"))

# Summarise means and sds of unscaled data
unscaled.means.sds <- as.data.frame(t(volcglacdata.inf.sub %>% summarise(mean(glacier.zmed, na.rm = TRUE), mean(glacier.slope, na.rm = TRUE),
                                                                         mean(glacier.thick.median, na.rm = TRUE), mean(climate.pca.2017.2018, na.rm = TRUE), 
                                                                         sd(glacier.zmed, na.rm = TRUE), sd(glacier.slope, na.rm = TRUE),
                                                                         sd(glacier.thick.median, na.rm = TRUE), sd(climate.pca.2017.2018, na.rm = TRUE))))
print(unscaled.means.sds)

# Scale and centre numerical variables
volcglacdata.inf.sub.scaled <- (as.data.frame(scale(volcglacdata.inf.sub)))

# Summarise means and sds of scaled data
scaled.means.sds <- as.data.frame(t(volcglacdata.inf.sub.scaled %>% summarise(mean(glacier.zmed, na.rm = TRUE), mean(glacier.slope, na.rm = TRUE),
                                                                              mean(glacier.thick.median, na.rm = TRUE), mean(climate.pca.2017.2018, na.rm = TRUE), 
                                                                              sd(glacier.zmed, na.rm = TRUE), sd(glacier.slope, na.rm = TRUE),
                                                                              sd(glacier.thick.median, na.rm = TRUE), sd(climate.pca.2017.2018, na.rm = TRUE))))
print(scaled.means.sds) # variable means are 0 and standard deviations are 1

# Add factor variables and unscaled dependent variable to scaled independent variables data frame
volcglacdata.inf.sub.scaled.add <- cbind(volcglacdata.inf.sub.scaled, volcglacdata.inf$glacier.type, volcglacdata.inf$glacier.vel.median)

# Rename added columns
colnames(volcglacdata.inf.sub.scaled.add)[which(names(volcglacdata.inf.sub.scaled.add) == "volcglacdata.inf$glacier.type")] <- "glacier.type"
colnames(volcglacdata.inf.sub.scaled.add)[which(names(volcglacdata.inf.sub.scaled.add) == "volcglacdata.inf$glacier.vel.median")] <- "glacier.vel.median"

# Drop data frame rows containing any NAs
volcglacdata.inf.sub.scaled.add.nadrop <- drop_na(volcglacdata.inf.sub.scaled.add)

# Parameterise 1st linear model (LM1)
LM1 <- lm(log10(glacier.vel.median) ~ glacier.type + glacier.thick.median + glacier.zmed + climate.pca.2017.2018,
          data = volcglacdata.inf.sub.scaled.add.nadrop)

# Output LM1 results
print(summary(LM1), digits = 4)

# Calculate LM1 confidence intervals
print(confint(LM1), digits = 2)

# Calculate LM1 AIC
AIC(LM1)

# Output LM1 result as table
LM1rep <- report(LM1)
as.report_table(LM1rep, summary = TRUE)

# Calculate LM1 stadardised residuals
LM1.sresid <- (LM1$residuals - mean(LM1$residuals))/sd(LM1$residuals)

# Assumption: standardised residuals should be normally distributed
hist(LM1.sresid, breaks = 100, freq = F)
lines(density(LM1.sresid, adjust = 1))
skewness(LM1.sresid)
qqPlot(LM1.sresid)

# Assumption: variances of standardised residuals should be homogenous (i.e. no heteroscedasity)
plot(LM1.sresid ~ LM1$fitted.values, pch = 20, cex = 1, cex.lab = 1)   # Plot shows some evidence of heteroscedasity in the dataset
plot(LM1.sresid ~ volcglacdata.inf.sub.scaled.add.nadrop$glacier.type)
plot(LM1.sresid ~ volcglacdata.inf.sub.scaled.add.nadrop$glacier.thick.median)   # Plot suggests source of heteroscedasity is the glacier thickness
plot(LM1.sresid ~ volcglacdata.inf.sub.scaled.add.nadrop$glacier.zmed)
plot(LM1.sresid ~ volcglacdata.inf.sub.scaled.add.nadrop$climate.pca.2017.2018)

# Assumption: the model is not biased by unduly influential observations
LM1.influence <- influence.measures(LM1)
View(LM1.influence$infmat)   # All cooks.d values < 1, therefore no unduly influential observations
LM1.CD <- cooks.distance(LM1)
LM1.CD.sresid <- resid(LM1, type = "pearson")
plot(LM1.CD~LM1.CD.sresid)   # NB. LM.CD is cooks.d

# Create new dataframe including log10 transformation of glacier thickness
volcglacdata.inf.sub.logthick <- cbind(volcglacdata.inf.sub, log10(volcglacdata.inf.sub$glacier.thick.median))
colnames(volcglacdata.inf.sub.logthick)[which(names(volcglacdata.inf.sub.logthick) == "log10(volcglacdata.inf.sub$glacier.thick.median)")] <- "glacier.thick.median.log10"

# Scale and centre numerical variables
volcglacdata.inf.sub.logthick.scaled <- (as.data.frame(scale(volcglacdata.inf.sub.logthick)))

# Summarise means and sds of scaled data
scaled.means.sds.logthick <- as.data.frame(t(volcglacdata.inf.sub.logthick.scaled %>% summarise(mean(glacier.zmed, na.rm = TRUE), mean(glacier.slope, na.rm = TRUE),
                                                                                                mean(glacier.thick.median.log10, na.rm = TRUE), mean(climate.pca.2017.2018, na.rm = TRUE), 
                                                                                                sd(glacier.zmed, na.rm = TRUE), sd(glacier.slope, na.rm = TRUE),
                                                                                                sd(glacier.thick.median.log10, na.rm = TRUE), sd(climate.pca.2017.2018, na.rm = TRUE))))
print(scaled.means.sds.logthick) # variable means are 0 and standard deviations are 1

# Add factor variables and unscaled dependent variable to scaled independent variables data frame
volcglacdata.inf.sub.logthick.scaled.add <- cbind(volcglacdata.inf.sub.logthick.scaled, volcglacdata.inf$glacier.type, volcglacdata.inf$glacier.vel.median)

# Rename added columns
colnames(volcglacdata.inf.sub.logthick.scaled.add)[which(names(volcglacdata.inf.sub.logthick.scaled.add) == "volcglacdata.inf$glacier.type")] <- "glacier.type"
colnames(volcglacdata.inf.sub.logthick.scaled.add)[which(names(volcglacdata.inf.sub.logthick.scaled.add) == "volcglacdata.inf$glacier.vel.median")] <- "glacier.vel.median"

# Drop data frame rows containing any NAs
volcglacdata.inf.sub.logthick.scaled.add.nadrop <- drop_na(volcglacdata.inf.sub.logthick.scaled.add)

# Parameterise 2nd linear model (LM2) with log10 glacier thickness
LM2 <- lm(log10(glacier.vel.median) ~ glacier.type + glacier.thick.median.log10 + glacier.zmed + climate.pca.2017.2018,
          data = volcglacdata.inf.sub.logthick.scaled.add.nadrop)

# Output LM2 results
print(summary(LM2), digits = 5)

# Calculate LM2 confidence intervals
print(confint(LM2), digits = 2)

# Calculate LM2 AIC
AIC(LM2)

# Output LM2 result as table
LM2rep <- report(LM2)
as.report_table(LM2rep, summary = TRUE)

# Calculate LM2 stadardised residuals
LM2.sresid <- (LM2$residuals - mean(LM2$residuals))/sd(LM2$residuals)

# Assumption: standardised residuals should be normally distributed
hist(LM2.sresid, breaks = 100, freq = F)
lines(density(LM2.sresid, adjust = 1))
skewness(LM2.sresid)
qqPlot(LM2.sresid)

# Assumption: variances of standardised residuals should be homogenous (i.e. no heteroscedasity)
plot(LM2.sresid ~ LM2$fitted.values, pch = 20, cex = 1, cex.lab = 1)   # Plot still shows some evidence of heteroscedasity in the dataset
plot(LM2.sresid ~ volcglacdata.inf.sub.logthick.scaled.add.nadrop$glacier.type)
plot(LM2.sresid ~ volcglacdata.inf.sub.logthick.scaled.add.nadrop$glacier.thick.median.log10)   # Plot suggests still some heteroscedasity in glacier thickness
plot(LM2.sresid ~ volcglacdata.inf.sub.logthick.scaled.add.nadrop$glacier.zmed)
plot(LM2.sresid ~ volcglacdata.inf.sub.logthick.scaled.add.nadrop$climate.pca.2017.2018)

# Assumption: the model is not biased by unduly influential observations
LM2.influence <- influence.measures(LM2)
View(LM2.influence$infmat)   # NB. Sort by cook.d column. All cooks.d values < 1, therefore no unduly influential observations
LM2.CD <- cooks.distance(LM2)
LM2.CD.sresid <- resid(LM2, type = "pearson")
plot(LM2.CD~LM2.CD.sresid)   # NB. LM.CD is cooks.d

# Generate Robust Linear Model (RLM) to investigate effect of outliers on model fit
LM2.robust<- rlm(log10(glacier.vel.median) ~ glacier.type + glacier.thick.median.log10 + glacier.zmed + climate.pca.2017.2018,
                 data = volcglacdata.inf.sub.logthick.scaled.add.nadrop)

# Output RLM results
print(summary(LM2.robust), digits = 3)

# Calculate RLM confidence intervals
print(confint.default(LM2.robust), digits = 2)

# Calculate RLM AIC
AIC(LM2.robust)

# Output RLM result as table
LM2.robust.rep <- report(LM2.robust)
as.report_table(LM2.robust.rep, summary = TRUE)   # NB. LM2 and RLM coefficients are consistent, so outliers have minimal effect on model fit

# Generate Generalised Least Squares Estimation Model (GLS) to address heteroscedasticity in LM2
GLS.baseline <- gls(log10(glacier.vel.median) ~ glacier.type + glacier.thick.median.log10 + glacier.zmed + climate.pca.2017.2018,
                    data = volcglacdata.inf.sub.logthick.scaled.add.nadrop)   # NB. Produces same coefficients as LM2 (because no weighting added yet)

GLS.power <- gls(log10(glacier.vel.median) ~ glacier.type + glacier.thick.median.log10 + glacier.zmed + climate.pca.2017.2018,
                 weights = varPower(form = ~glacier.thick.median.log10),   # Includes a power weighting
                 data = volcglacdata.inf.sub.logthick.scaled.add.nadrop)

GLS.exp <- gls(log10(glacier.vel.median) ~ glacier.type + glacier.thick.median.log10 + glacier.zmed + climate.pca.2017.2018,
               weights = varExp(form = ~glacier.thick.median.log10),   # Includes an exponential weighting
               data = volcglacdata.inf.sub.logthick.scaled.add.nadrop)

GLS.constpower <- gls(log10(glacier.vel.median) ~ glacier.type + glacier.thick.median.log10 + glacier.zmed + climate.pca.2017.2018,
                      weights = varConstPower(form = ~glacier.thick.median.log10),   # Includes a constant and power weighting
                      data = volcglacdata.inf.sub.logthick.scaled.add.nadrop)

# Compare GLS AICs
AIC(GLS.baseline, GLS.power, GLS.exp, GLS.constpower)

# Calculate GLS stadardised residuals
GLS.baseline.sresid <- (GLS.baseline$residuals - mean(GLS.baseline$residuals))/sd(GLS.baseline$residuals)
GLS.power.sresid <- (GLS.power$residuals - mean(GLS.power$residuals))/sd(GLS.power$residuals)
GLS.exp.sresid <- (GLS.exp$residuals - mean(GLS.exp$residuals))/sd(GLS.exp$residuals)
GLS.constpower.sresid <- (GLS.constpower$residuals - mean(GLS.constpower$residuals))/sd(GLS.constpower$residuals)

# Compare GLS standardised residual skewness
skewness(GLS.baseline.sresid)
skewness(GLS.power.sresid)
skewness(GLS.exp.sresid)
skewness(GLS.constpower.sresid)

# Plot standarised GLS residuals against glacier thickness to visualise effect on homoscedasticity
plot(GLS.baseline.sresid ~ volcglacdata.inf.sub.logthick.scaled.add.nadrop$glacier.thick.median.log10)
plot(GLS.power.sresid ~ volcglacdata.inf.sub.logthick.scaled.add.nadrop$glacier.thick.median.log10)
plot(GLS.exp.sresid ~ volcglacdata.inf.sub.logthick.scaled.add.nadrop$glacier.thick.median.log10)
plot(GLS.constpower.sresid ~ volcglacdata.inf.sub.logthick.scaled.add.nadrop$glacier.thick.median.log10)
# NB. No GLM weightings are able to substantially address heteroscedasticity in the glacier thickness residuals, so use outputs of LM2





#### Volcanic Controls on Glacier Velocity: Descriptive Statistics ####

# Create new dataframe with that only includes volcanic glaciers with a median velocity value
onlyvolc <- volcglacdata.inf %>% filter(glacier.type == "volcanic glacier")

# Tables, violin plots, ANOVAs and Tukey tests of volcanic glacier velocity by volcanic property
# Volcano type
as.data.frame(onlyvolc %>% group_by(volcano.type) %>% summarise(mean(glacier.vel.median), median(glacier.vel.median), sd(glacier.vel.median), n())) %>% mutate_if(is.numeric, ~round(., 2))
as.data.frame(onlyvolc %>% group_by(volcano.type) %>% count(volcano.number) %>% summarise (n()))   # Counts no. of UNIQUE volcanoes by volcano property
violin.volc.type <- ggplot(onlyvolc, aes(x = volcano.type, y = glacier.vel.median, fill = volcano.type)) +
  geom_violin() +
  scale_y_continuous(trans = "log10", labels = scales::comma, minor_breaks = minor_breaks_log(10)) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "Volcano type", y = expression(Median~velocity~~(ma^-1))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("Stratovolcano","Shield", "Complex", "Caldera", "Pyroclastic cone", "Subglacial", "Volcanic field")) +
  guides(fill="none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, vjust = 0), axis.title.x = element_text(size = 12, vjust = 0)) +
  scale_fill_brewer(palette = "Pastel1", direction=-1)
violin.volc.type
anova.volcano.type <- lm(glacier.vel.median ~ volcano.type, data = onlyvolc)   # NB. type 2 ANOVA for unbalanced groups
Anova(anova.volcano.type, type = 2)   ### p = 0.01, thus ANOVA suggests a significant difference between two or more group means...
aov.volcano.type <- aov(glacier.vel.median ~ volcano.type, data = onlyvolc)
TukeyHSD(aov.volcano.type)   ### ...but no significant differences when using Tukey HSD to compare all possible pairs of means

# Tectonic setting
as.data.frame(onlyvolc %>% group_by(tectonic.setting) %>% summarise(mean(glacier.vel.median), median(glacier.vel.median), sd(glacier.vel.median), n())) %>% mutate_if(is.numeric, ~round(., 2))
as.data.frame(onlyvolc %>% group_by(tectonic.setting) %>% count(volcano.number) %>% summarise (n()))   # Counts no. of UNIQUE volcanoes by volcano property
violin.tectonic.setting <- ggplot(onlyvolc, aes(x=tectonic.setting, y=glacier.vel.median, fill=tectonic.setting)) +
  geom_violin() +
  scale_y_continuous(trans = "log10", labels = scales::comma, minor_breaks = minor_breaks_log(10)) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "Tectonic setting", y = expression(Median~velocity~~(ma^-1))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("Subduct. zone (cont. crust)","Subduct. zone (intermed. crust)", "Rift zone (ocean crust)", "Intraplate (cont. crust)", "Intraplate (ocean crust)")) +
  guides(fill="none") +
  theme_minimal() +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, vjust = 0), axis.title.x = element_text(size = 12, vjust = 0)) +
  scale_fill_brewer(palette = "Pastel1", direction=-1)
violin.tectonic.setting
anova.tectonic.setting <- lm(glacier.vel.median ~ tectonic.setting, data = onlyvolc)   # NB. type 2 ANOVA for unbalanced groups
Anova(anova.tectonic.setting, type = 2)   ### p = <0.001, thus ANOVA suggests a significant difference between two or more group means...
aov.tectonic.setting <- aov(glacier.vel.median ~ tectonic.setting, data = onlyvolc)
TukeyHSD(aov.tectonic.setting)

# Dominant rock type
as.data.frame(onlyvolc %>% group_by(dominant.rock.type) %>% summarise(mean(glacier.vel.median), median(glacier.vel.median), sd(glacier.vel.median), n())) %>% mutate_if(is.numeric, ~round(., 2))
as.data.frame(onlyvolc %>% group_by(dominant.rock.type) %>% count(volcano.number) %>% summarise (n()))   # Counts no. of UNIQUE volcanoes by volcano property
violin.dominant.rock.type <- ggplot(onlyvolc, aes(x=dominant.rock.type, y=glacier.vel.median, fill=dominant.rock.type)) +
  geom_violin() +
  scale_y_continuous(trans = "log10", labels = scales::comma, minor_breaks = minor_breaks_log(10)) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  geom_boxplot(width = 0.2, colour = "grey25", alpha = 0.1) +
  labs(x = "Dominant rock type", y = expression(Median~velocity~~(ma^-1))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  guides(fill="none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, vjust = 0), axis.title.x = element_text(size = 12, vjust = 0)) +
  scale_fill_brewer(palette = "Pastel1", direction=-1)
violin.dominant.rock.type
anova.dominant.rock.type <- lm(glacier.vel.median ~ dominant.rock.type, data = onlyvolc)   # NB. type 2 ANOVA for unbalanced groups
Anova(anova.dominant.rock.type, type = 2)   ### p = <0.001, thus ANOVA suggests a significant difference between two or more group means...
aov.dominant.rock.type <- aov(glacier.vel.median ~ dominant.rock.type, data = onlyvolc)
TukeyHSD(aov.dominant.rock.type)

# Volcano region
as.data.frame(onlyvolc %>% group_by(volcano.region) %>% summarise(mean(glacier.vel.median), median(glacier.vel.median), sd(glacier.vel.median), n())) %>% mutate_if(is.numeric, ~round(., 2))
as.data.frame(onlyvolc %>% group_by(volcano.region) %>% count(volcano.number) %>% summarise (n()))   # Counts no. of UNIQUE volcanoes by volcano property
violin.volcano.region <- ggplot(onlyvolc, aes(x=volcano.region, y=glacier.vel.median, fill=volcano.region)) +
  geom_violin() +
  scale_y_continuous(trans = "log10", labels = scales::comma, minor_breaks = minor_breaks_log(10)) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "Volcano region", y = expression(Median~velocity~~(ma^-1))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("Alaska", "Canada & West. USA", "Iceland & Arctic Ocean", "Kamchatka & Mainland Asia", "Mediterranean & West. Asia", "Mid. East & Indian Ocean", "South America")) +
  guides(fill="none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, vjust = 0), axis.title.x = element_text(size = 12, vjust = 0)) +
  scale_fill_brewer(palette = "Pastel1", direction=-1)
violin.volcano.region
anova.volcano.region <- lm(glacier.vel.median ~ volcano.region, data = onlyvolc)   # NB. type 2 ANOVA for unbalanced groups
Anova(anova.volcano.region, type = 2)   ### p = <0.001, thus ANOVA suggests a significant difference between two or more group means...
aov.volcano.region <- aov(glacier.vel.median ~ volcano.region, data = onlyvolc)
TukeyHSD(aov.volcano.region)

# Glacier distance to volcano
as.data.frame(onlyvolc %>% group_by(glacier.distance.to.vent) %>% summarise(mean(glacier.vel.median), median(glacier.vel.median), sd(glacier.vel.median), n())) %>% mutate_if(is.numeric, ~round(., 2))
as.data.frame(onlyvolc %>% group_by(glacier.distance.to.vent) %>% count(volcano.number) %>% summarise (n()))   # Counts no. of UNIQUE volcanoes by volcano property
violin.glacier.distance.to.vent <- ggplot(onlyvolc, aes(x=glacier.distance.to.vent, y=glacier.vel.median, fill=glacier.distance.to.vent)) +
  geom_violin(alpha = 0.9) +
  scale_y_continuous(trans = "log10", labels = scales::comma, minor_breaks = minor_breaks_log(10)) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "Glacier proximity to volcano (km)", y = expression(Median~velocity~~(ma^-1))) +
  scale_x_discrete(labels=c("≤ 1.0", "> 1.0 - 2.5", "> 2.5 - 5.0")) +
  guides(fill="none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "grey",  linewidth = 1), panel.grid.major.x = element_line(color = "grey",  linewidth = 0.1),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "solid", linewidth = 0.1),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11, vjust = 0), axis.title.x = element_text(size = 11, vjust = 0)) +
  scale_fill_brewer(palette = "OrRd", direction=-1)
violin.glacier.distance.to.vent
anova.glacier.distance.to.vent <- lm(glacier.vel.median ~ glacier.distance.to.vent, data = onlyvolc)   # NB. type 2 ANOVA for unbalanced groups
Anova(anova.glacier.distance.to.vent, type = 2)   ### p = 0.002, thus ANOVA suggests a significant difference between two or more group means...
aov.glacier.distance.to.vent <- aov(glacier.vel.median ~ glacier.distance.to.vent, data = onlyvolc)
TukeyHSD(aov.glacier.distance.to.vent)

# No. of volcanoes within 5 km of glacier
as.data.frame(onlyvolc %>% group_by(glacier.no.volc.centres) %>% summarise(mean(glacier.vel.median), median(glacier.vel.median), sd(glacier.vel.median), n())) %>% mutate_if(is.numeric, ~round(., 2))
as.data.frame(onlyvolc %>% group_by(glacier.no.volc.centres) %>% count(volcano.number) %>% summarise (n()))   # Counts no. of UNIQUE volcanoes by volcano property
onlyvolc$glacier.no.volc.centres <- as.factor(onlyvolc$glacier.no.volc.centres)   # NB. change variable from numeric to factor
violin.glacier.no.volc.centres <- ggplot(onlyvolc, aes(x=glacier.no.volc.centres, y=glacier.vel.median, fill=glacier.no.volc.centres)) +
  geom_violin() +
  scale_y_continuous(trans = "log10", labels = scales::comma, minor_breaks = minor_breaks_log(10)) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "Number of volcanoes within a 5 km radius of glacier", y = expression(Median~velocity~~(ma^-1))) +
  scale_x_discrete() +
  guides(fill="none") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title.y=element_text(size=12, vjust=2), axis.title.x=element_text(size=12, vjust=-1)) +
  scale_fill_brewer(palette = "Pastel1", direction=-1)
violin.glacier.no.volc.centres
anova.glacier.no.volc.centres <- lm(glacier.vel.median ~ glacier.no.volc.centres, data = onlyvolc)   # NB. type 2 ANOVA for unbalanced groups
Anova(anova.glacier.no.volc.centres, type = 2)   ### p = 0.175, thus ANOVA suggests no significant differences group means

# Eruption evidence
as.data.frame(onlyvolc %>% group_by(eruption.evidence) %>% summarise(mean(glacier.vel.median), median(glacier.vel.median), sd(glacier.vel.median), n())) %>% mutate_if(is.numeric, ~round(., 2))
as.data.frame(onlyvolc %>% group_by(eruption.evidence) %>% count(volcano.number) %>% summarise (n()))   # Counts no. of UNIQUE volcanoes by volcano property
violin.eruption.evidence <- ggplot(onlyvolc, aes(x=eruption.evidence, y=glacier.vel.median, fill=eruption.evidence)) +
  geom_violin() +
  scale_y_continuous(trans = "log10", labels = scales::comma, minor_breaks = minor_breaks_log(10)) +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE) +
  geom_boxplot(width = 0.3, colour = "grey25", alpha = 0.1) +
  labs(x = "Eruption evidence", y = expression(Median~velocity~~(ma^-1))) +
  scale_x_discrete() +
  guides(fill="none") +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, vjust = 0), axis.title.x = element_text(size = 12, vjust = 0)) +
  scale_fill_brewer(palette = "Pastel1", direction=-1)
violin.eruption.evidence
anova.eruption.evidence <- lm(glacier.vel.median ~ eruption.evidence, data = onlyvolc)   # NB. type 2 ANOVA for unbalanced groups
Anova(anova.eruption.evidence, type = 2)   ### p = <0.001, thus ANOVA suggests a significant difference between two or more group means...
aov.eruption.evidence <- aov(glacier.vel.median ~ eruption.evidence, data = onlyvolc)
TukeyHSD(aov.eruption.evidence)

# Last eruption (aBP; BP = 2019) 
onlyvolc %>% filter(eruption.abp != "NA") %>% summarise(n())
scatter.last.eruption <- ggplot(data = onlyvolc, mapping = aes(x = eruption.abp, y = glacier.vel.median)) +
  scale_x_continuous(trans=compose_trans("log10", "reverse")) +
  scale_y_log10(limits = c(0.5,500)) +
  geom_pointdensity() +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  scale_color_viridis() +
  labs(x = "Years since last eruption (years before 2019)", y = expression(Median~velocity~~(ma^-1))) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, vjust = 0), axis.title.x = element_text(size = 12, vjust = 0))
as.data.frame(onlyvolc %>% group_by(eruption.abp) %>% count(volcano.number) %>% summarise (n()))   # Counts no. of UNIQUE volcanoes by volcano property
scatter.last.eruption
lm.last.eruption <- lm(glacier.vel.median ~ eruption.abp, data = onlyvolc)
summary(lm.last.eruption)

# No. Holocene eruptions
onlyvolc %>% filter(glacier.vel.median != "NA") %>% summarise(n())
scatter.no.holocene.eruptions <- ggplot(data = onlyvolc, mapping = aes(x = no.holocene.eruptions, y = glacier.vel.median)) +
  scale_x_continuous() +
  scale_y_log10(limits = c(0.5, 500)) +
  geom_pointdensity() +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  scale_color_viridis() +
  labs(x = "Number of Holocene eruptions", y = expression(Median~velocity~~(ma^-1))) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dashed"),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, vjust = 0), axis.title.x = element_text(size = 12, vjust = 0))
as.data.frame(onlyvolc %>% group_by(no.holocene.eruptions) %>% count(volcano.number) %>% summarise (n()))   # Counts no. of UNIQUE volcanoes by volcano property
scatter.no.holocene.eruptions
lm.no.holocene.eruptions <- lm(glacier.vel.median ~ no.holocene.eruptions, data = onlyvolc)
summary(lm.no.holocene.eruptions)





#### Volcanic Controls on Glacier Velocity: Inferential Statistics ####

# Histogram of dependent variable to determine distribution
hist(onlyvolc$glacier.vel.median, breaks = 100)   # I.e. dependent variable possesses an inverse gaussian distribution
skewness(onlyvolc$glacier.vel.median)
qqPlot(onlyvolc$glacier.vel.median)

# Apply log10 transformation to correct for inverse gaussian distribution
hist(log10(onlyvolc$glacier.vel.median), breaks = 100)
skewness(log10(onlyvolc$glacier.vel.median))
qqPlot(log10(onlyvolc$glacier.vel.median))

# Check for sparse variables in volcglacdata.inf dataframe
nearZeroVar(onlyvolc)
colnames(onlyvolc)

# Create new dataframe to check for multicollinearity in potential model variables
onlyvolc.cor <- onlyvolc

# Convert factor variables to numeric variables for correlation matrix
onlyvolc.cor$glacier.no.volc.centres <- as.numeric(as.character(onlyvolc.cor$glacier.no.volc.centres))

# Correlation matrix for all potential numeric model variables:
correlation.volc.variables <- cor(onlyvolc.cor[c("glacier.thick.median", "glacier.zmed", "climate.pca.2017.2018",
                                                 "glacier.no.volc.centres", "eruption.abp", "no.holocene.eruptions")],
                                  use = "pairwise.complete.obs")
corrplot(cor(correlation.volc.variables), method = "number", type = "full", col = COL2('RdBu', 10), tl.srt = 45, cl.ratio = 0.2, tl.cex = 0.8, number.cex = 0.8)

# Count n for potential model variables
sum(is.na(onlyvolc$glacier.vel.median))
sum(is.na(onlyvolc$glacier.thick.median))
sum(is.na(onlyvolc$glacier.zmed))
sum(is.na(onlyvolc$climate.pca.2017.2018))
sum(is.na(onlyvolc$volcano.type))
sum(is.na(onlyvolc$tectonic.setting))
sum(is.na(onlyvolc$dominant.rock.type))
sum(is.na(onlyvolc$volcano.region))
sum(is.na(onlyvolc$glacier.distance.to.vent))
sum(is.na(onlyvolc$glacier.no.volc.centres))
sum(is.na(onlyvolc$eruption.evidence))
sum(is.na(onlyvolc$eruption.abp))
sum(is.na(onlyvolc$no.holocene.eruptions))
sum(is.na(onlyvolc$volcano.number))

# Correlation matrix for likely potential numeric model variables:
correlation.volc.variables.likely <- cor(onlyvolc.cor[c("glacier.thick.median", "glacier.zmed", "climate.pca.2017.2018",
                                                        "glacier.no.volc.centres", "no.holocene.eruptions")],
                                         use = "pairwise.complete.obs")
LMM.corrplot.likely <- cor(correlation.volc.variables.likely)
colnames(LMM.corrplot.likely) <- c("Median\n thickness", "Median elevation", "Climate PC", "No. volcanoes within 5 km", "No. Holocene eruptions")
rownames(LMM.corrplot.likely) <- c("Median\n thickness", "Median elevation", "Climate PC", "No. volcanoes within 5 km", "No. Holocene eruptions")
corrplot(LMM.corrplot.likely, method = "color", diag = FALSE, outline = TRUE, type = "lower", tl.col = "black", addCoef.col = TRUE, col = COL2('RdBu', 100), tl.srt = 90, cl.ratio = 0.2, tl.cex = 0.8, number.cex = 0.8)

# Variance Inflation Factor (VIF) to check for multicollinearity between all potential numerical and factor variables
onlyvolc.vif <- lm(glacier.vel.median ~ glacier.thick.median + glacier.zmed + climate.pca.2017.2018 +
                     volcano.type + tectonic.setting + dominant.rock.type + volcano.region + glacier.distance.to.vent +
                     glacier.no.volc.centres + eruption.evidence + eruption.abp + no.holocene.eruptions,  data = onlyvolc)
as.data.frame(vif(onlyvolc.vif))   # NB. Returns error: "there are aliased coefficients in the model" - i.e. correlation between some independent variables

# Run alias to identify correlated variables
alias(lm(glacier.vel.median ~ glacier.thick.median + glacier.zmed + climate.pca.2017.2018 +
           volcano.type + tectonic.setting + dominant.rock.type + volcano.region + glacier.distance.to.vent +
           glacier.no.volc.centres + eruption.evidence + eruption.abp + no.holocene.eruptions,  data = onlyvolc))   ### Tectonic setting and volcano region are correlated, only retain the former in the model

# Re-run VIF without volcano region variable
onlyvolc.vif2 <- lm(glacier.vel.median ~ glacier.thick.median + glacier.zmed + climate.pca.2017.2018 +
                      volcano.type + tectonic.setting + dominant.rock.type + glacier.distance.to.vent +
                      glacier.no.volc.centres + eruption.evidence + eruption.abp + no.holocene.eruptions,  data = onlyvolc)
as.data.frame(vif(onlyvolc.vif2))   # All values <5, so no multicollinearity between factor variables and numerical variables (use GVIF for numerical variables and GVIF^(1/2*Df) for categorical variables)

# Variance Inflation Factor (VIF) for check for multicollinearity between numerical and factor variables in final LMM
onlyvolc.vif.final <- lm(glacier.vel.median ~ glacier.thick.median + glacier.zmed + climate.pca.2017.2018 +
                           volcano.type + tectonic.setting + glacier.distance.to.vent + glacier.no.volc.centres + no.holocene.eruptions,  data = onlyvolc)
as.data.frame(vif(onlyvolc.vif.final))   # All values <5, so no multicollinearity between factor variables and numerical variables (use GVIF for numerical variables and GVIF^(1/2*Df) for categorical variables)

# Subset dataframe with likely model variables that require scaling and centring
onlyvolc.sub <- subset(onlyvolc, select = c("glacier.thick.median", "glacier.zmed", "climate.pca.2017.2018",
                                            "eruption.abp", "no.holocene.eruptions"))

# Add log10 of glacier thickness in case required to address normality of the model residuals ############ Can take out if not required
onlyvolc.sub$glacier.thick.median.log <- log10(onlyvolc.sub$glacier.thick.median)

# Summarise means and sds of unscaled data
unscaled.volconly.means.sds <- as.data.frame(t(onlyvolc.sub %>% summarise(mean(glacier.thick.median, na.rm = T), mean(glacier.thick.median.log, na.rm = T),
                                                                          mean(glacier.zmed, na.rm = T), mean(climate.pca.2017.2018, na.rm = T),
                                                                          mean(eruption.abp, na.rm = T), mean(no.holocene.eruptions, na.rm = T),
                                                                          sd(glacier.thick.median, na.rm = T), sd(glacier.thick.median.log, na.rm = T),
                                                                          sd(glacier.zmed, na.rm = T), sd(climate.pca.2017.2018, na.rm = T),
                                                                          sd(eruption.abp, na.rm = T), sd(no.holocene.eruptions, na.rm = T))))
print(unscaled.volconly.means.sds)

# Scale and centre numerical variables
onlyvolc.sub.scaled <- (as.data.frame(scale(onlyvolc.sub)))

# Summarise means and sds of scaled data
scaled.volconly.means.sds <- as.data.frame(t(onlyvolc.sub.scaled %>% summarise(mean(glacier.thick.median, na.rm = T), mean(glacier.thick.median.log, na.rm = T),
                                                                               mean(glacier.zmed, na.rm = T), mean(climate.pca.2017.2018, na.rm = T),
                                                                               mean(eruption.abp, na.rm = T), mean(no.holocene.eruptions, na.rm = T),
                                                                               sd(glacier.thick.median, na.rm = T), sd(glacier.thick.median.log, na.rm = T),
                                                                               sd(glacier.zmed, na.rm = T), sd(climate.pca.2017.2018, na.rm = T),
                                                                               sd(eruption.abp, na.rm = T), sd(no.holocene.eruptions, na.rm = T))))
print(scaled.volconly.means.sds)   # variable means are now 0 and standard deviations are 1

# Add factor variables, random effect (volcano number), and unscaled dependent variable to scaled independent variable data frame
onlyvolc.sub.scaled.add <- cbind(onlyvolc.sub.scaled, onlyvolc$volcano.type, onlyvolc$tectonic.setting, onlyvolc$dominant.rock.type,
                                 onlyvolc$volcano.region, onlyvolc$glacier.distance.to.vent, onlyvolc$glacier.no.volc.centres,
                                 onlyvolc$eruption.evidence, onlyvolc$volcano.number, onlyvolc$glacier.vel.median, onlyvolc$rgi.id)

# Rename added columns
colnames(onlyvolc.sub.scaled.add)[which(names(onlyvolc.sub.scaled.add) == "onlyvolc$volcano.type")] <- "volcano.type"
colnames(onlyvolc.sub.scaled.add)[which(names(onlyvolc.sub.scaled.add) == "onlyvolc$tectonic.setting")] <- "tectonic.setting"
colnames(onlyvolc.sub.scaled.add)[which(names(onlyvolc.sub.scaled.add) == "onlyvolc$dominant.rock.type")] <- "dominant.rock.type"
colnames(onlyvolc.sub.scaled.add)[which(names(onlyvolc.sub.scaled.add) == "onlyvolc$volcano.region")] <- "volcano.region"
colnames(onlyvolc.sub.scaled.add)[which(names(onlyvolc.sub.scaled.add) == "onlyvolc$glacier.distance.to.vent")] <- "glacier.distance.to.vent"
colnames(onlyvolc.sub.scaled.add)[which(names(onlyvolc.sub.scaled.add) == "onlyvolc$glacier.no.volc.centres")] <- "glacier.no.volc.centres"
colnames(onlyvolc.sub.scaled.add)[which(names(onlyvolc.sub.scaled.add) == "onlyvolc$eruption.evidence")] <- "eruption.evidence"
colnames(onlyvolc.sub.scaled.add)[which(names(onlyvolc.sub.scaled.add) == "onlyvolc$volcano.number")] <- "volcano.number"
colnames(onlyvolc.sub.scaled.add)[which(names(onlyvolc.sub.scaled.add) == "onlyvolc$glacier.vel.median")] <- "glacier.vel.median"
colnames(onlyvolc.sub.scaled.add)[which(names(onlyvolc.sub.scaled.add) == "onlyvolc$rgi.id")] <- "rgi.id"

# Parameterise 1st linear mixed model (LMM1) and test assumptions
onlyvolc.sub.scaled.add.lmm1 <- onlyvolc.sub.scaled.add
lmm1 <- lmer(log10(glacier.vel.median) ~ glacier.thick.median + glacier.zmed + climate.pca.2017.2018 + volcano.type +
               tectonic.setting + glacier.no.volc.centres + glacier.distance.to.vent + no.holocene.eruptions + (1 | volcano.number),
             REML = T, data = onlyvolc.sub.scaled.add.lmm1)
summary(lmm1)

# Calculate lmm1 confidence intervals
print(confint(lmm1), digits = 2)

# lmm1 anova
anova(lmm1)

# Calculate lmm1 AIC
AIC(lmm1)
AICc(lmm1)

# Output lmm1 result as table
lmm1.rep <- report(lmm1)
as.report_table(lmm1.rep, summary = TRUE)

# Calculate lmm1 stadardised residuals
lmm1.sresid <- (resid(lmm1) - mean(resid(lmm1)))/sd(resid(lmm1))

# Assumption: standardised residuals should be normally distributed
hist(lmm1.sresid, breaks = 100, freq = F)
lines(density(lmm1.sresid, adjust = 1))
skewness(lmm1.sresid)
qqPlot(lmm1.sresid)

# Assumption: variances of standardised residuals should be homogenous (i.e. no heteroscedasity)
plot(lmm1.sresid ~ fitted.values(lmm1), pch = 20, cex = 1, cex.lab = 1)   # Plot suggests slight heteroscedasticity in the dataset
onlyvolc.sub.scaled.add.forlmm1plots <-  within(onlyvolc.sub.scaled.add.lmm1, rm("eruption.abp", "eruption.evidence",  
                                                                                 "dominant.rock.type"))
onlyvolc.sub.scaled.add.forlmm1plots <-  na.omit(onlyvolc.sub.scaled.add.forlmm1plots)
plot(lmm1.sresid ~ onlyvolc.sub.scaled.add.forlmm1plots$glacier.thick.median)   # Median thickness is source of heteroscedasticity
plot(lmm1.sresid ~ onlyvolc.sub.scaled.add.forlmm1plots$glacier.zmed)
plot(lmm1.sresid ~ onlyvolc.sub.scaled.add.forlmm1plots$climate.pca.2017.2018)
plot(lmm1.sresid ~ onlyvolc.sub.scaled.add.forlmm1plots$tectonic.setting)
plot(lmm1.sresid ~ onlyvolc.sub.scaled.add.forlmm1plots$glacier.distance.to.vent)
plot(lmm1.sresid ~ onlyvolc.sub.scaled.add.forlmm1plots$glacier.no.volc.centres)
plot(lmm1.sresid ~ onlyvolc.sub.scaled.add.forlmm1plots$no.holocene.eruptions)

# Assumption: the dataset does not contain serial auto-correlation
durbinWatsonTest(resid(lmm1))   # p > 0.05 therefore no auto-correlation

# Assumption: the model is not biased by unduly influential observations
lmm1.CD <- cooks.distance(lmm1)
lmm1.CD.sresid <- resid(lmm1, type = "pearson")
plot(lmm1.CD~lmm1.CD.sresid)   # NB. One datapoint CD > 1 and therefore unduly influential
lmm1.CDtable <- as.data.frame(cooks.distance(lmm1))
View(lmm1.CDtable)   # NB. Datapoint #6970 (glacier ID = RGI60-01.08367) is unduly influential

# Parameterise 2nd linear mixed model (LMM2) and test assumptions (model omits glacier "RGI60-01.08367")
onlyvolc.sub.scaled.add.lmm2 <- onlyvolc.sub.scaled.add
onlyvolc.sub.scaled.add.lmm2$glacier.thick.median[which(onlyvolc.sub.scaled.add.lmm2$rgi.id == "RGI60-01.08367")] <- NA   # NB. removes glacier RGI60-01.08367 from dataset
lmm2 <- lmer(log10(glacier.vel.median) ~ glacier.thick.median + glacier.zmed + climate.pca.2017.2018 + volcano.type +
                tectonic.setting + glacier.no.volc.centres + glacier.distance.to.vent + no.holocene.eruptions + (1 | volcano.number),
              REML = T, data = onlyvolc.sub.scaled.add.lmm2)
summary(lmm2)

# Calculate lmm2 confidence intervals
print(confint(lmm2), digits = 2)

# lmm2 anova
anova(lmm2)

# Calculate lmm2 AIC
AIC(lmm2)
AICc(lmm2)

# Output lmm2 result as table
lmm2.rep <- report(lmm2)
as.report_table(lmm2.rep, summary = TRUE)

# Calculate lmm2 stadardised residuals
lmm2.sresid <- (resid(lmm2) - mean(resid(lmm2)))/sd(resid(lmm2))

# Assumption: standardised residuals should be normally distributed
hist(lmm2.sresid, breaks = 100, freq = F)
lines(density(lmm2.sresid, adjust = 1))
skewness(lmm2.sresid)
qqPlot(lmm2.sresid)

# Assumption: variances of standardised residuals should be homogenous (i.e. no heteroscedasity)
plot(lmm2.sresid ~ fitted.values(lmm2), pch = 20, cex = 1, cex.lab = 1)   # Plot suggests slight heteroscedasticity in the dataset
onlyvolc.sub.scaled.add.lmm2$glacier.thick.median[which(onlyvolc.sub.scaled.add.lmm2$rgi.id == "RGI60-01.08367")] <- NA   # NB. removes RGI60-01.08367 from dataset
onlyvolc.sub.scaled.add.forlmm2plots <-  within(onlyvolc.sub.scaled.add.lmm2, rm("eruption.abp", "eruption.evidence",  
                                                                                   "dominant.rock.type"))
onlyvolc.sub.scaled.add.forlmm2plots <-  na.omit(onlyvolc.sub.scaled.add.forlmm2plots)
plot(lmm2.sresid ~ onlyvolc.sub.scaled.add.forlmm2plots$glacier.thick.median)   # Median thickness is source of heteroscedasticity
plot(lmm2.sresid ~ onlyvolc.sub.scaled.add.forlmm2plots$glacier.zmed)
plot(lmm2.sresid ~ onlyvolc.sub.scaled.add.forlmm2plots$climate.pca.2017.2018)
plot(lmm2.sresid ~ onlyvolc.sub.scaled.add.forlmm2plots$tectonic.setting)
plot(lmm2.sresid ~ onlyvolc.sub.scaled.add.forlmm2plots$glacier.distance.to.vent)
plot(lmm2.sresid ~ onlyvolc.sub.scaled.add.forlmm2plots$glacier.no.volc.centres)
plot(lmm2.sresid ~ onlyvolc.sub.scaled.add.forlmm2plots$no.holocene.eruptions)

# Assumption: the dataset does not contain serial auto-correlation
durbinWatsonTest(resid(lmm2))   # p > 0.05 therefore no auto-correlation

# Assumption: the model is not biased by unduly influential observations
lmm2.CD <- cooks.distance(lmm2)
lmm2.CD.sresid <- resid(lmm2, type = "pearson")
plot(lmm2.CD~lmm2.CD.sresid)   # All cooks.d values (LM.CD) < 1, therefore no unduly influential observations
lmm2.CDtable <- as.data.frame(cooks.distance(lmm2))
View(lmm2.CDtable)

# LMM comparison
# Compare LMM R2s   # NB. R2m = variance explained by fixed factors only, R2c = variance explained by fixed factors AND random factor (volcano ID)
r.squaredGLMM(lmm1)
r.squaredGLMM(lmm2)

# Compare LMM AICs
AIC(lmm1, lmm2)

# Compare LMM RSEs
summary(lmm1)$sigma
summary(lmm2)$sigma

# Generate Robust Linear Mixed Model (RLMM) to investigate effect of outliers on model fit
lmm2.robust <- rlmer(log10(glacier.vel.median) ~ glacier.thick.median + glacier.zmed + climate.pca.2017.2018 +
                        volcano.type + tectonic.setting + glacier.no.volc.centres + glacier.distance.to.vent +
                        no.holocene.eruptions + (1 | volcano.number),
                      data = onlyvolc.sub.scaled.add.lmm2)

# Output RLMM results
print(summary(lmm2.robust), digits = 6)   # NB. LMM2b and RLMM (lmm2b.robust) coefficients consistent, so outliers have minimal effect on model fit

# Calculate RLMM confidence intervals
confint.rlmerMod <- function(lmm2.robust, level = 0.95) {
  beta <- fixef(lmm2.robust)
  parm <- names(beta)
  se <- sqrt(diag(vcov(lmm2.robust)))
  z <- qnorm((1 + level) / 2)
  ctab <- cbind(beta - (z * se), 
                beta + (z * se))
  colnames(ctab) <- c(paste(100 * ((1 - level) / 2), '%'),
                      paste(100 * ((1 + level) / 2), '%'))
  return(ctab[parm, ])}
confint(lmm2.robust)

# Calculate RLMM RSE (NB. no p values, R2 or AIC available for RLMMs)
summary(lmm2.robust)$sigma