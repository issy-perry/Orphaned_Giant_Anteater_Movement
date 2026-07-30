# load packages
library(ctmm)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(terra) # main package we need to use for dowloading, cropping, and saving our rasters
library(sf) # getting the extents of datasets
library(weights)
library(tidyverse)

# load tel data and UDs
load("./DATA/Wild_raised_tel_data/Data_telemetry_RR.rda") #wild-raised 
load("./DATA/Orphaned_tel_data/Data_telemetry_RR.rda") #orphaned 
load("./RESULTS/AKDEs/UDs_orphan_RR.rda") # orphaned 
load("./RESULTS/AKDEs/UDs_wild_RR.rda") # wild-raised

# separate tel data and UDs by years monitored
# wild-raised
DATA_17 <- DATA_wild[c("Anthony", "Bumpus", "Cate", "Christoffer", "Elaine", "Jackson", 
                       "Kyle", "Little_Rick", "Makao", "Puji")]
DATA_18 <- DATA_wild[c("Alexander", "Annie", "Beto", "Hannah", "Jane", "Larry", "Luigi",
                       "Margaret", "Maria", "Reid", "Rodolfo", "Sheron", "Thomas")] 
# orphaned
DATA_19 <- DATA_orphan[c("Arya", "Dumbo_1", "Dumbo_2")]
DATA_21 <- DATA_orphan[c("Tim_2")]
DATA_22 <- DATA_orphan[c("Colete", "Heather", "Juju_2", "Mulan", "Peter", "Tim_3")] 
DATA_23 <- DATA_orphan[c("Bella", "George", "Nancy")] 
DATA_24 <- DATA_orphan[c("Dom")] 
#wild-raised
AKDE_17 <- AKDE_wild[c("Anthony", "Bumpus", "Cate", "Christoffer", "Elaine", "Jackson", 
                       "Kyle", "Little_Rick", "Makao", "Puji")] 
AKDE_18 <- AKDE_wild[c("Alexander", "Annie", "Beto", "Hannah", "Jane", "Larry", "Luigi",
                       "Margaret", "Maria", "Reid", "Rodolfo", "Sheron", "Thomas")] 
# orphaned
AKDE_19 <- AKDE_orphan[c("Arya", "Dumbo_1", "Dumbo_2")] 
AKDE_21 <- AKDE_orphan[c("Tim_2")] 
AKDE_22 <- AKDE_orphan[c("Colete", "Heather", "Juju_2", "Mulan", "Peter", "Tim_3")] 
AKDE_23 <- AKDE_orphan[c("Bella", "George", "Nancy")] 
AKDE_24 <- AKDE_orphan[c("Dom")] 

# call raster downloading script to run 
source("./SCRIPTS/000 Downloading Land Coverage Rasters.R")

# 2017 ----
USE_17 <- data.frame()
for(i in 1:length(DATA_17)){
  # subset individual
  AKDES <- AKDE_17[[i]]
  # rename raster
  land_types <- cover_2017
  
  HR <- rast(raster(AKDES, DF = "PMF"))
  HR2 <- project(HR, crs(land_types), res = res(land_types))
  HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
  HR <- project(HR, crs(land_types), res = res(land_types))
  HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df$layer <- HR.df$layer/sum(HR.df$layer)
  #Extract habitat values
  HR.df$land_class <- raster::extract(land_types, HR.df[,1:2])[,2]
  HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29", "9")] <- "Forest"
  HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
  HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
  HR.df$land_class[HR.df$land_class %in% c("24","25","30", "75")] <- "Development"
  HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
  
  #Use the home range PDF to calculate the weighted proportions of time spent the different land class types
  PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
  PROPS2 <- data.frame(class = names(PROPS),
                       proportion = as.numeric(PROPS))
  PROPS <- data.frame(t(PROPS2))[2,]
  names(PROPS) <- PROPS2$class
  #make dataframe
  USE <- data.frame(binomial = "Myrmecophaga_tridactyla")
  USE$ID <- AKDES@info$identity
  #bind results to dataframe
  USE <- cbind(USE,PROPS)
  # add to total dataframe
  USE_17 <- bind_rows(USE_17, USE)
} # end of loop

# 2018 ----
USE_18 <- data.frame()
for(i in 1:length(DATA_18)){
  # subset individual
  AKDES <- AKDE_18[[i]]
  # rename raster
  land_types <- cover_2018
  
  HR <- rast(raster(AKDES, DF = "PMF"))
  HR2 <- project(HR, crs(land_types), res = res(land_types))
  HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
  HR <- project(HR, crs(land_types), res = res(land_types))
  HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df$layer <- HR.df$layer/sum(HR.df$layer)
  #Extract habitat values
  HR.df$land_class <- raster::extract(land_types, HR.df[,1:2])[,2]
  HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29", "9")] <- "Forest"
  HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
  HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
  HR.df$land_class[HR.df$land_class %in% c("24","25","30", "75")] <- "Development"
  HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
  
  #Use the home range PDF to calculate the weighted proportions of time spent the different land class types
  PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
  PROPS2 <- data.frame(class = names(PROPS),
                       proportion = as.numeric(PROPS))
  PROPS <- data.frame(t(PROPS2))[2,]
  names(PROPS) <- PROPS2$class
  #make dataframe
  USE <- data.frame(binomial = "Myrmecophaga_tridactyla")
  USE$ID <- AKDES@info$identity
  #bind results to dataframe
  USE <- cbind(USE,PROPS)
  # add to total dataframe
  USE_18 <- bind_rows(USE_18, USE)
} # end of loop


# 2019 ----
USE_19 <- data.frame()
for(i in 1:length(DATA_19)){
  # subset individual
  AKDES <- AKDE_19[[i]]
  # rename raster
  land_types <- cover_2019
  
  HR <- rast(raster(AKDES, DF = "PMF"))
  HR2 <- project(HR, crs(land_types), res = res(land_types))
  HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
  HR <- project(HR, crs(land_types), res = res(land_types))
  HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df$layer <- HR.df$layer/sum(HR.df$layer)
  #Extract habitat values
  HR.df$land_class <- raster::extract(land_types, HR.df[,1:2])[,2]
  HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29", "9")] <- "Forest"
  HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
  HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
  HR.df$land_class[HR.df$land_class %in% c("24","25","30", "75")] <- "Development"
  HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
  
  #Use the home range PDF to calculate the weighted proportions of time spent the different land class types
  PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
  PROPS2 <- data.frame(class = names(PROPS),
                       proportion = as.numeric(PROPS))
  PROPS <- data.frame(t(PROPS2))[2,]
  names(PROPS) <- PROPS2$class
  #make dataframe
  USE <- data.frame(binomial = "Myrmecophaga_tridactyla")
  USE$ID <- AKDES@info$identity
  #bind results to dataframe
  USE <- cbind(USE,PROPS)
  # add to total dataframe
  USE_19 <- bind_rows(USE_19, USE)
} # end of loop


# 2021 ----
USE_21 <- data.frame()
for(i in 1:length(DATA_21)){
  # subset individual
  AKDES <- AKDE_21[[i]]
  # rename raster
  land_types <- cover_2021
  
  HR <- rast(raster(AKDES, DF = "PMF"))
  HR2 <- project(HR, crs(land_types), res = res(land_types))
  HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
  HR <- project(HR, crs(land_types), res = res(land_types))
  HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df$layer <- HR.df$layer/sum(HR.df$layer)
  #Extract habitat values
  HR.df$land_class <- raster::extract(land_types, HR.df[,1:2])[,2]
  HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29", "9")] <- "Forest"
  HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
  HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
  HR.df$land_class[HR.df$land_class %in% c("24","25","30", "75")] <- "Development"
  HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
  
  #Use the home range PDF to calculate the weighted proportions of time spent the different land class types
  PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
  PROPS2 <- data.frame(class = names(PROPS),
                       proportion = as.numeric(PROPS))
  PROPS <- data.frame(t(PROPS2))[2,]
  names(PROPS) <- PROPS2$class
  #make dataframe
  USE <- data.frame(binomial = "Myrmecophaga_tridactyla")
  USE$ID <- AKDES@info$identity
  #bind results to dataframe
  USE <- cbind(USE,PROPS)
  # add to total dataframe
  USE_21 <- bind_rows(USE_21, USE)
} # end of loop


# 2022 ----
USE_22 <- data.frame()
for(i in 1:length(DATA_22)){
  # subset individual
  AKDES <- AKDE_22[[i]]
  # rename raster
  land_types <- cover_2022
  
  HR <- rast(raster(AKDES, DF = "PMF"))
  HR2 <- project(HR, crs(land_types), res = res(land_types))
  HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
  HR <- project(HR, crs(land_types), res = res(land_types))
  HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df$layer <- HR.df$layer/sum(HR.df$layer)
  #Extract habitat values
  HR.df$land_class <- raster::extract(land_types, HR.df[,1:2])[,2]
  HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29", "9")] <- "Forest"
  HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
  HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
  HR.df$land_class[HR.df$land_class %in% c("24","25","30", "75")] <- "Development"
  HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
  
  #Use the home range PDF to calculate the weighted proportions of time spent the different land class types
  PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
  PROPS2 <- data.frame(class = names(PROPS),
                       proportion = as.numeric(PROPS))
  PROPS <- data.frame(t(PROPS2))[2,]
  names(PROPS) <- PROPS2$class
  #make dataframe
  USE <- data.frame(binomial = "Myrmecophaga_tridactyla")
  USE$ID <- AKDES@info$identity
  #bind results to dataframe
  USE <- cbind(USE,PROPS)
  # add to total dataframe
  USE_22 <- bind_rows(USE_22, USE)
} # end of loop


# 2023 ----
USE_23 <- data.frame()
for(i in 1:length(DATA_23)){
  # subset individual
  AKDES <- AKDE_23[[i]]
  # rename raster
  land_types <- cover_2023
  
  HR <- rast(raster(AKDES, DF = "PMF"))
  HR2 <- project(HR, crs(land_types), res = res(land_types))
  HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
  HR <- project(HR, crs(land_types), res = res(land_types))
  HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df$layer <- HR.df$layer/sum(HR.df$layer)
  #Extract habitat values
  HR.df$land_class <- raster::extract(land_types, HR.df[,1:2])[,2]
  HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29", "9")] <- "Forest"
  HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
  HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
  HR.df$land_class[HR.df$land_class %in% c("24","25","30", "75")] <- "Development"
  HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
  
  #Use the home range PDF to calculate the weighted proportions of time spent the different land class types
  PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
  PROPS2 <- data.frame(class = names(PROPS),
                       proportion = as.numeric(PROPS))
  PROPS <- data.frame(t(PROPS2))[2,]
  names(PROPS) <- PROPS2$class
  #make dataframe
  USE <- data.frame(binomial = "Myrmecophaga_tridactyla")
  USE$ID <- AKDES@info$identity
  #bind results to dataframe
  USE <- cbind(USE,PROPS)
  # add to total dataframe
  USE_23 <- bind_rows(USE_23, USE)
} # end of loop


# 2024 ----
USE_24 <- data.frame()
for(i in 1:length(DATA_24)){
  # subset individual
  AKDES <- AKDE_24[[i]]
  # rename raster
  land_types <- cover_2024
  
  HR <- rast(raster(AKDES, DF = "PMF"))
  HR2 <- project(HR, crs(land_types), res = res(land_types))
  HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
  HR <- project(HR, crs(land_types), res = res(land_types))
  HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
  #Renormalize
  HR.df$layer <- HR.df$layer/sum(HR.df$layer)
  #Extract habitat values
  HR.df$land_class <- raster::extract(land_types, HR.df[,1:2])[,2]
  HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29", "9")] <- "Forest"
  HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
  HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
  HR.df$land_class[HR.df$land_class %in% c("24","25","30", "75")] <- "Development"
  HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
  
  #Use the home range PDF to calculate the weighted proportions of time spent the different land class types
  PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
  PROPS2 <- data.frame(class = names(PROPS),
                       proportion = as.numeric(PROPS))
  PROPS <- data.frame(t(PROPS2))[2,]
  names(PROPS) <- PROPS2$class
  #make dataframe
  USE <- data.frame(binomial = "Myrmecophaga_tridactyla")
  USE$ID <- AKDES@info$identity
  #bind results to dataframe
  USE <- cbind(USE,PROPS)
  # add to total dataframe
  USE_24 <- bind_rows(USE_24, USE)
} # end of loop


# bind dataframes together
Wild_use <- bind_rows(USE_17, USE_18)
Orphan_use <- bind_rows(USE_19, USE_21, USE_22, USE_23, USE_24)

# replace NAs with 0s
Wild_use[is.na(Wild_use)] <- 0
Orphan_use[is.na(Orphan_use)] <- 0

# add information about population to each
Wild_use$Site <- "Wild-raised"
Orphan_use$Site <- "Orphaned"

# ensure that columns are numerical
Wild_use[c(3:7)] <- lapply(Wild_use[c(3:7)], as.numeric)
Orphan_use[c(3:7)] <- lapply(Orphan_use[c(3:7)], as.numeric)

# means of proportions across years wild-raised site
wild_mean <- data.frame()
wild_mean[1, "Agriculture"] <- mean(Wild_use[["Agriculture"]])
wild_mean[1, "Development"] <- mean(Wild_use[["Development"]])
wild_mean[1, "Forest"] <- mean(Wild_use[["Forest"]])
wild_mean[1, "Pasture"] <- mean(Wild_use[["Pasture"]])
wild_mean[1, "Water"] <- mean(Wild_use[["Water"]])
orphan_mean <- data.frame()
orphan_mean[1, "Agriculture"] <- mean(Orphan_use[["Agriculture"]])
orphan_mean[1, "Development"] <- mean(Orphan_use[["Development"]])
orphan_mean[1, "Forest"] <- mean(Orphan_use[["Forest"]])
orphan_mean[1, "Pasture"] <- mean(Orphan_use[["Pasture"]])
orphan_mean[1, "Water"] <- mean(Orphan_use[["Water"]])

# add information about population to each
wild_mean$Site <- "Wild-raised"
orphan_mean$Site <- "Orphaned"
total_use <- rbind(wild_mean, orphan_mean)

#pivot to make plotting easier
total_use_piv <- total_use %>% pivot_longer(cols=c("Water", "Agriculture", "Development", "Pasture", "Forest"),
                                            names_to = "Land_Type",
                                            values_to = "Percentage")

# plot Use
total_use_piv$Site <- factor(total_use_piv$Site, levels = c("Wild-raised", "Orphaned"))


use_plot <- ggplot(total_use_piv, aes(x = Land_Type, y = Percentage, fill = Site)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual("Population", 
                    values = c("Wild-raised" = "#F5AD5B", "Orphaned" = "#79DECD"),
                    name = "Population") +
  labs(title = "", x = "", y = "Percentage of Use") +
  
  theme(plot.title = element_text(size = 12, hjust = 0.065),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA), 
        legend.position = "right", 
        legend.title = element_text(size=12, family = "sans", face = "bold", hjust = 0.5),
        legend.text = element_text(size = 11, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        #strip.text = element_text(face = "bold", size = 12, family = "sans"),
        axis.line = element_line(color = "darkgray", linewidth = 0.5),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 12)) 

plot(use_plot)

ggsave("./FIGURES/Revised/Supplementary_Figure_use.png", width = 8.5, height = 5.5, units = "in", 
       dpi = 300, bg = "transparent")

