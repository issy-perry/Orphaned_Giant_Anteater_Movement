# 6 Window Analysis
# This script includes the following:
#     1. Extract land types from rasters to create covariates for RSFs
#     2. Fit RSFs using tel data and AKDEs against each covariate

#The moving window code is primarily adapted from Kat Chhen's analysis for her masters thesis 
#This code is based on Ryan Gill's moving window (which in turn is based on Dr. Michael Noonan's code)



# load packages 
library(ctmm) # includes functions for movement models, UDs, and RSFs
library(terra) # extracting land types from Spatrasters prior to converting to raster classes
library(raster) # convert Spatrasters to raster classes for RSFs
library(lutz) #working with times
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(lubridate) #working with dates
library(tictoc) #seeing how long each individual takes to run
library(sf)
library(crayon) #adding colored bands to printed announcements for tictoc
library(weights) # For calculating the weighted percentage of the land class types

# load data
load("./DATA/Wild_raised_tel_data/Data_telemetry_RR.rda") #wild-raised PICKLES
load("./DATA/Orphaned_tel_data/Data_telemetry_RR.rda") #orphaned PICKLES

# set time periods for window
# window was determined by META output of tau position for range-resident orphans (2.1 days est) and wild-raised (1.4 days est) 
# the average value between the two was 1.75. 1.75 multiplied by 6 equaled roughly 10 days. This is sufficient to reduce temporal autocorrelation by over 95%.
dt <- 2 %#% "day" # window shifts 2 days
win <- 10 %#% "day" # window covers 10 days (based on average tau position value extracted using ctmm::meta() (1.669) multiplied by 6)

#create folders to hold results
folder_list <- c("Fits", 
                 "AKDEs",
                 "Mean_Speed",
                 "Land_Use")


#set directory path
dir_path <- "./RESULTS/Window/" 

# create every folder in the folder list in the directory path
for (folder in folder_list) {
  dir.create(paste0(dir_path, folder), 
             recursive = TRUE, showWarnings = TRUE)
}



# Separate Individuals based on year monitored ----
DATA_17 <- DATA_wild[c("Anthony", "Bumpus", "Cate", "Christoffer", "Elaine", "Jackson", "Kyle", "Little_Rick",
                       "Makao", "Puji")] 
DATA_18 <- DATA_wild[c("Alexander", "Annie", "Beto", "Hannah", "Jane", "Larry", "Luigi", "Margaret", "Maria",
                       "Reid", "Rodolfo", "Sheron", "Thomas")]  
# orphaned
DATA_19 <- DATA_orphan[c("Arya", "Dumbo_1", "Dumbo_2")]
#DATA_20 <- DATA_orphan[c("Tim_1")] 
DATA_21 <- DATA_orphan[c("Tim_2")]
DATA_22 <- DATA_orphan[c("Colete", "Heather", "Juju_2", "Mulan", "Peter", "Rita", "Tim_3")] 
DATA_23 <- DATA_orphan[c("Bella", "George", "Nancy")]
DATA_24 <- DATA_orphan[c("Dom")]

# call raster downloading script to run 
source("./SCRIPTS/000 Downloading Land Coverage Rasters.R")

# Window 2017 ----
# print which window year is being worked on
cat(bgGreen(paste0("window 2017")), "\n")

# for loop for window analysis
for(i in 1:length(DATA_17)){
  tic("window analysis 2017")
  cat(bgGreen(paste0("on iteration ", i)), "\n")
  # subset out an individual
  DATA <- DATA_17[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 2 days forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), 
                      origin = "1970-01-01", 
                      tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                  DATA$longitude[1], 
                                                  method = "accurate")) 
  
  # set up list to store
  wAKDEs <- list()
  FITS <- list()
  speed_mean <- list()
  #RSFs <- list()
  use_land <- list()
  
  
  
  # Analysis on the window segment ......................................
  for (j in 1:length(times)) {
    # extract data within the window segment
    SUBSET <- DATA[times[j] <= DATA$t & DATA$t <= times[j] + win,] # +win means window size (2 days)
    
    if (nrow(SUBSET) == 0) {
      cat("No data found for window section in iteration", j, "- moving on to the next iteration.\n")
      next
    }
    
    # get subset window start and end based on the recorded collar data
    WINDOW_START <- as.POSIXct(min(SUBSET$t), origin = "1970-01-01", 
                               tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    WINDOW_END <- as.POSIXct(max(SUBSET$t), origin = "1970-01-01", 
                             tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    
    # Indicate the iteration and window segment 
    cat(bgMagenta(paste((j), "of", length(times), "iterations. Window segment:",
                        WINDOW_START, "to", WINDOW_END,
                        "for anteater:", DATA@info[1]), "\n"))
    cat(paste0("Number of fixes in window segment subset: ", nrow(SUBSET), "\n"))
    
    # Process the subset if data is present
    tryCatch({
      cat(bgBlue(paste0("Movement Model  ")), "\n")
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 1, cores = 4))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        cat(bgMagenta(paste0("UD  ")), "\n")
        AKDES <- akde(SUBSET, FIT, weights = TRUE, cores = 4)
        
        #mean speeds
        cat(bgCyan(paste0("Speed  ")), "\n")
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = 4, trace = FALSE) #in m/s
        
        
        #weighted proportion of land use
        cat(bgGreen(paste0("Proportion of Use  ")), "\n")
        HR <- rast(raster(AKDES, DF = "PMF"))
        HR2 <- project(HR, crs(cover_2017), res = res(cover_2017))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(cover_2017), res = res(cover_2017))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        HR.df$land_class <- raster::extract(cover_2017, HR.df[,1:2])[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("9","1","3", "4","5","6","49","29")] <- "Forested"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        
        #Use the home range PDF to calculate the weighted percentage of time spent the different land class types
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
        
        # store results in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        #RSFs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- RSF
        use_land[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- USE
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rda for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  #save(RSFs, file = paste0(dir_path, "RSF_Fit/RSF_", DATA@info[1], ".rda"))
  save(use_land, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  
  # clean up environment
  rm(FIT, AKDES, SPEED, USE) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}
# clear up environment space
rm(FITS, wAKDEs, speed_mean, use_land, DATA_17)
gc()


# Window 2018 ----
# print which window year is being worked on
cat(bgGreen(paste0("window 2018")), "\n")

# for loop for window analysis
for(i in 1:length(DATA_18)){
  tic("window analysis 2018")
  cat(bgGreen(paste0("on iteration ", i)), "\n")
  # subset out an individual
  DATA <- DATA_18[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 2 days forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), 
                      origin = "1970-01-01", 
                      tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                  DATA$longitude[1], 
                                                  method = "accurate")) 
  
  # set up list to store
  wAKDEs <- list()
  FITS <- list()
  speed_mean <- list()
  #RSFs <- list()
  use_land <- list()
  
  # Analysis on the window segment ......................................
  for (j in 1:length(times)) {
    # extract data within the window segment
    SUBSET <- DATA[times[j] <= DATA$t & DATA$t <= times[j] + win,] # +win means window size (2 days)
    
    if (nrow(SUBSET) == 0) {
      cat("No data found for window section in iteration", j, "- moving on to the next iteration.\n")
      next
    }
    
    # get subset window start and end based on the recorded collar data
    WINDOW_START <- as.POSIXct(min(SUBSET$t), origin = "1970-01-01", 
                               tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    WINDOW_END <- as.POSIXct(max(SUBSET$t), origin = "1970-01-01", 
                             tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    
    # Indicate the iteration and window segment 
    cat(bgMagenta(paste((j), "of", length(times), "iterations. Window segment:",
                        WINDOW_START, "to", WINDOW_END,
                        "for anteater:", DATA@info[1]), "\n"))
    cat(paste0("Number of fixes in window segment subset: ", nrow(SUBSET), "\n"))
    
    # Process the subset if data is present
    tryCatch({
      cat(bgBlue(paste0("Movement Model  ")), "\n")
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 1, cores = 4))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        cat(bgMagenta(paste0("UD  ")), "\n")
        AKDES <- akde(SUBSET, FIT, weights = TRUE, cores = 4)
        
        #mean speeds
        cat(bgCyan(paste0("Speed  ")), "\n")
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = 4, trace = FALSE) #in m/s
        
        
        #weighted proportion of land use
        cat(bgGreen(paste0("Proportion of Use  ")), "\n")
        HR <- rast(raster(AKDES, DF = "PMF"))
        HR2 <- project(HR, crs(cover_2018), res = res(cover_2018))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(cover_2018), res = res(cover_2018))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        HR.df$land_class <- raster::extract(cover_2018, HR.df[,1:2])[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("9","1","3", "4","5","6","49","29")] <- "Forested"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        
        #Use the home range PDF to calculate the weighted percentage of time spent the different land class types
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
        
        # store results in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        #RSFs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- RSF
        use_land[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- USE
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rda for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  #save(RSFs, file = paste0(dir_path, "RSF_Fit/RSF_", DATA@info[1], ".rda"))
  save(use_land, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  
  # clean up environment
  rm(FIT, AKDES, SPEED, USE) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

# clear environment space
rm(FITS, wAKDEs, speed_mean, use_land, DATA_18)
gc()

# Window 2019 ----

# download raster
land_cover_2019 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2019.tif")
# crop raster based on extent of wild-raised population
land_types <- crop(land_cover_2019, eo)

# remove previous large raster 
rm(land_cover_2019)
gc()
# print which window year is being worked on
cat(bgGreen(paste0("window 2019")), "\n")

# for loop for window analysis
for(i in 1:length(DATA_19)){
  tic("window analysis 2019")
  cat(bgGreen(paste0("on iteration ", i)), "\n")
  # subset out an individual
  DATA <- DATA_19[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 2 days forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), 
                      origin = "1970-01-01", 
                      tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                  DATA$longitude[1], 
                                                  method = "accurate")) 
  
  # set up list to store
  wAKDEs <- list()
  FITS <- list()
  speed_mean <- list()
  #RSFs <- list()
  use_land <- list()
  
  # Analysis on the window segment ......................................
  for (j in 1:length(times)) {
    # extract data within the window segment
    SUBSET <- DATA[times[j] <= DATA$t & DATA$t <= times[j] + win,] # +win means window size (2 days)
    
    if (nrow(SUBSET) == 0) {
      cat("No data found for window section in iteration", j, "- moving on to the next iteration.\n")
      next
    }
    
    # get subset window start and end based on the recorded collar data
    WINDOW_START <- as.POSIXct(min(SUBSET$t), origin = "1970-01-01", 
                               tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    WINDOW_END <- as.POSIXct(max(SUBSET$t), origin = "1970-01-01", 
                             tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    
    # Indicate the iteration and window segment 
    cat(bgMagenta(paste((j), "of", length(times), "iterations. Window segment:",
                        WINDOW_START, "to", WINDOW_END,
                        "for anteater:", DATA@info[1]), "\n"))
    cat(paste0("Number of fixes in window segment subset: ", nrow(SUBSET), "\n"))
    
    # Process the subset if data is present
    tryCatch({
      cat(bgBlue(paste0("Movement Model  ")), "\n")
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 1, cores = 4))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        cat(bgMagenta(paste0("UD  ")), "\n")
        AKDES <- akde(SUBSET, FIT, weights = TRUE, cores = 4)
        
        #mean speeds
        cat(bgCyan(paste0("Speed  ")), "\n")
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = 4, trace = FALSE) #in m/s
        
        #weighted proportion of land use
        cat(bgGreen(paste0("Proportion of Use  ")), "\n")
        HR <- rast(raster(AKDES, DF = "PMF"))
        HR2 <- project(HR, crs(cover_2019), res = res(cover_2019))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(cover_2019), res = res(cover_2019))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        HR.df$land_class <- raster::extract(cover_2019, HR.df[,1:2])[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("9","1","3", "4","5","6","49","29")] <- "Forested"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        
        #Use the home range PDF to calculate the weighted percentage of time spent the different land class types
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
        
        # store results in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        #RSFs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- RSF
        use_land[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- USE
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rda for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  #save(RSFs, file = paste0(dir_path, "RSF_Fit/RSF_", DATA@info[1], ".rda"))
  save(use_land, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  
  # clean up environment
  rm(FIT, AKDES, SPEED, USE) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}


rm(FITS, wAKDEs, speed_mean, use_land, DATA_19)



# Window 2021 ----

# download raster
land_cover_2021 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2021.tif")
# crop raster based on extent of wild-raised population
land_types <- crop(land_cover_2021, eo)

# remove previous large raster 
rm(land_cover_2021)
gc()
# print which window year is being worked on
cat(bgGreen(paste0("window 2021")), "\n")

# for loop for window analysis
for(i in 1:length(DATA_21)){
  tic("window analysis 2021")
  cat(bgGreen(paste0("on iteration ", i)), "\n")
  # subset out an individual
  DATA <- DATA_21[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 2 days forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), 
                      origin = "1970-01-01", 
                      tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                  DATA$longitude[1], 
                                                  method = "accurate")) 
  
  # set up list to store
  wAKDEs <- list()
  FITS <- list()
  speed_mean <- list()
  #RSFs <- list()
  use_land <- list()
  
  # Analysis on the window segment ......................................
  for (j in 1:length(times)) {
    # extract data within the window segment
    SUBSET <- DATA[times[j] <= DATA$t & DATA$t <= times[j] + win,] # +win means window size (2 days)
    
    if (nrow(SUBSET) == 0) {
      cat("No data found for window section in iteration", j, "- moving on to the next iteration.\n")
      next
    }
    
    # get subset window start and end based on the recorded collar data
    WINDOW_START <- as.POSIXct(min(SUBSET$t), origin = "1970-01-01", 
                               tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    WINDOW_END <- as.POSIXct(max(SUBSET$t), origin = "1970-01-01", 
                             tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    
    # Indicate the iteration and window segment 
    cat(bgMagenta(paste((j), "of", length(times), "iterations. Window segment:",
                        WINDOW_START, "to", WINDOW_END,
                        "for anteater:", DATA@info[1]), "\n"))
    cat(paste0("Number of fixes in window segment subset: ", nrow(SUBSET), "\n"))
    
    # Process the subset if data is present
    tryCatch({
      cat(bgBlue(paste0("Movement Model  ")), "\n")
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 1, cores = 4))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        cat(bgMagenta(paste0("UD  ")), "\n")
        AKDES <- akde(SUBSET, FIT, weights = TRUE, cores = 4)
        
        #mean speeds
        cat(bgCyan(paste0("Speed  ")), "\n")
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = 4, trace = FALSE) #in m/s
        
        
        #weighted proportion of land use
        cat(bgGreen(paste0("Proportion of Use  ")), "\n")
        HR <- rast(raster(AKDES, DF = "PMF"))
        HR2 <- project(HR, crs(cover_2021), res = res(cover_2021))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(cover_2021), res = res(cover_2021))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        HR.df$land_class <- raster::extract(cover_2021, HR.df[,1:2])[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("9","1","3", "4","5","6","49","29")] <- "Forested"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        
        #Use the home range PDF to calculate the weighted percentage of time spent the different land class types
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
        
        # store results in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        #RSFs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- RSF
        use_land[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- USE
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rda for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  #save(RSFs, file = paste0(dir_path, "RSF_Fit/RSF_", DATA@info[1], ".rda"))
  save(use_land, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  
  # clean up environment
  rm(FIT, AKDES, SPEED, USE) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

# clear environment space
rm(FITS, wAKDEs, speed_mean, use_land, DATA_21)
gc()


# Window 2022 ----

# download raster
land_cover_2022 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2022.tif")
# crop raster based on extent of wild-raised population
land_types <- crop(land_cover_2022, eo)

# remove previous large raster 
rm(land_cover_2022)
gc()
# print which window year is being worked on
cat(bgGreen(paste0("window 2022")), "\n")

# for loop for window analysis
for(i in 1:length(DATA_22)){
  tic("window analysis 2022")
  cat(bgGreen(paste0("on iteration ", i)), "\n")
  # subset out an individual
  DATA <- DATA_22[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 2 days forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), 
                      origin = "1970-01-01", 
                      tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                  DATA$longitude[1], 
                                                  method = "accurate")) 
  
  # set up list to store
  wAKDEs <- list()
  FITS <- list()
  speed_mean <- list()
  #RSFs <- list()
  use_land <- list()
  
  # Analysis on the window segment ......................................
  for (j in 1:length(times)) {
    # extract data within the window segment
    SUBSET <- DATA[times[j] <= DATA$t & DATA$t <= times[j] + win,] # +win means window size (2 days)
    
    if (nrow(SUBSET) == 0) {
      cat("No data found for window section in iteration", j, "- moving on to the next iteration.\n")
      next
    }
    
    # get subset window start and end based on the recorded collar data
    WINDOW_START <- as.POSIXct(min(SUBSET$t), origin = "1970-01-01", 
                               tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    WINDOW_END <- as.POSIXct(max(SUBSET$t), origin = "1970-01-01", 
                             tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    
    # Indicate the iteration and window segment 
    cat(bgMagenta(paste((j), "of", length(times), "iterations. Window segment:",
                        WINDOW_START, "to", WINDOW_END,
                        "for anteater:", DATA@info[1]), "\n"))
    cat(paste0("Number of fixes in window segment subset: ", nrow(SUBSET), "\n"))
    
    # Process the subset if data is present
    tryCatch({
      cat(bgBlue(paste0("Movement Model  ")), "\n")
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 1, cores = 4))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        cat(bgMagenta(paste0("UD  ")), "\n")
        AKDES <- akde(SUBSET, FIT, weights = TRUE, cores = 4)
        
        #mean speeds
        cat(bgCyan(paste0("Speed  ")), "\n")
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = 4, trace = FALSE) #in m/s
        
        
        #weighted proportion of land use
        cat(bgGreen(paste0("Proportion of Use  ")), "\n")
        HR <- rast(raster(AKDES, DF = "PMF"))
        HR2 <- project(HR, crs(cover_2022), res = res(cover_2022))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(cover_2022), res = res(cover_2022))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        HR.df$land_class <- raster::extract(cover_2022, HR.df[,1:2])[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("9","1","3", "4","5","6","49","29")] <- "Forested"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        
        #Use the home range PDF to calculate the weighted percentage of time spent the different land class types
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
        
        # store results in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        #RSFs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- RSF
        use_land[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- USE
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rda for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  #save(RSFs, file = paste0(dir_path, "RSF_Fit/RSF_", DATA@info[1], ".rda"))
  save(use_land, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  
  # clean up environment
  rm(FIT, AKDES, SPEED, USE) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

# clear environment space
rm(FITS, wAKDEs, speed_mean, use_land, DATA_22)
gc()



# Window 2023 ----

# download raster
land_cover_2023 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2023.tif")
# crop raster based on extent of wild-raised population
land_types <- crop(land_cover_2023, eo)

# remove previous large raster 
rm(land_cover_2023)
gc()
# print which window year is being worked on
cat(bgGreen(paste0("window 2023")), "\n")

# for loop for window analysis
for(i in 1:length(DATA_23)){
  tic("window analysis 2023")
  cat(bgGreen(paste0("on iteration ", i)), "\n")
  # subset out an individual
  DATA <- DATA_23[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 2 days forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), 
                      origin = "1970-01-01", 
                      tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                  DATA$longitude[1], 
                                                  method = "accurate")) 
  
  # set up list to store
  wAKDEs <- list()
  FITS <- list()
  speed_mean <- list()
  #RSFs <- list()
  use_land <- list()
  
  
  # Analysis on the window segment ......................................
  for (j in 1:length(times)) {
    # extract data within the window segment
    SUBSET <- DATA[times[j] <= DATA$t & DATA$t <= times[j] + win,] # +win means window size (2 days)
    
    if (nrow(SUBSET) == 0) {
      cat("No data found for window section in iteration", j, "- moving on to the next iteration.\n")
      next
    }
    
    # get subset window start and end based on the recorded collar data
    WINDOW_START <- as.POSIXct(min(SUBSET$t), origin = "1970-01-01", 
                               tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    WINDOW_END <- as.POSIXct(max(SUBSET$t), origin = "1970-01-01", 
                             tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    
    # Indicate the iteration and window segment 
    cat(bgMagenta(paste((j), "of", length(times), "iterations. Window segment:",
                        WINDOW_START, "to", WINDOW_END,
                        "for anteater:", DATA@info[1]), "\n"))
    cat(paste0("Number of fixes in window segment subset: ", nrow(SUBSET), "\n"))
    
    # Process the subset if data is present
    tryCatch({
      cat(bgBlue(paste0("Movement Model  ")), "\n")
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 1, cores = 4))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        cat(bgMagenta(paste0("UD  ")), "\n")
        AKDES <- akde(SUBSET, FIT, weights = TRUE, cores = 4)
        
        #mean speeds
        cat(bgCyan(paste0("Speed  ")), "\n")
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = 4, trace = FALSE) #in m/s
        
        
        #weighted proportion of land use
        cat(bgGreen(paste0("Proportion of Use  ")), "\n")
        HR <- rast(raster(AKDES, DF = "PMF"))
        HR2 <- project(HR, crs(cover_2023), res = res(cover_2023))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(cover_2023), res = res(cover_2023))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        HR.df$land_class <- raster::extract(cover_2023, HR.df[,1:2])[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("9","1","3", "4","5","6","49","29")] <- "Forested"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        
        #Use the home range PDF to calculate the weighted percentage of time spent the different land class types
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
    
        
        # store results in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        #RSFs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- RSF
        use_land[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- USE
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rda for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  #save(RSFs, file = paste0(dir_path, "RSF_Fit/RSF_", DATA@info[1], ".rda"))
  save(use_land, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  
  # clean up environment
  rm(FIT, AKDES, SPEED, USE) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

# clear environment space
rm(FITS, wAKDEs, speed_mean, use_land, DATA_23)
gc()



# Window 2024 ----

# download raster
land_cover_2024 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2024.tif")
# crop raster based on extent of wild-raised population
land_types <- crop(land_cover_2024, eo)

# remove previous large raster 
rm(land_cover_2024)
gc()
# print which window year is being worked on
cat(bgGreen(paste0("window 2024")), "\n")

# for loop for window analysis
for(i in 1:length(DATA_24)){
  tic("window analysis 2024")
  cat(bgGreen(paste0("on iteration ", i)), "\n")
  # subset out an individual
  DATA <- DATA_24[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 2 days forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), 
                      origin = "1970-01-01", 
                      tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                  DATA$longitude[1], 
                                                  method = "accurate")) 
  
  # set up list to store
  wAKDEs <- list()
  FITS <- list()
  speed_mean <- list()
  #RSFs <- list()
  use_land <- list()
  
  # Analysis on the window segment ......................................
  for (j in 1:length(times)) {
    # extract data within the window segment
    SUBSET <- DATA[times[j] <= DATA$t & DATA$t <= times[j] + win,] # +win means window size (2 days)
    
    if (nrow(SUBSET) == 0) {
      cat("No data found for window section in iteration", j, "- moving on to the next iteration.\n")
      next
    }
    
    # get subset window start and end based on the recorded collar data
    WINDOW_START <- as.POSIXct(min(SUBSET$t), origin = "1970-01-01", 
                               tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    WINDOW_END <- as.POSIXct(max(SUBSET$t), origin = "1970-01-01", 
                             tz = lutz::tz_lookup_coords(SUBSET$latitude[1], SUBSET$longitude[1], method = "fast")) 
    
    # Indicate the iteration and window segment 
    cat(bgMagenta(paste((j), "of", length(times), "iterations. Window segment:",
                        WINDOW_START, "to", WINDOW_END,
                        "for anteater:", DATA@info[1]), "\n"))
    cat(paste0("Number of fixes in window segment subset: ", nrow(SUBSET), "\n"))
    
    # Process the subset if data is present
    tryCatch({
      cat(bgBlue(paste0("Movement Model  ")), "\n")
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 1, cores = 4))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        cat(bgMagenta(paste0("UD  ")), "\n")
        AKDES <- akde(SUBSET, FIT, weights = TRUE, cores = 4)
        
        #mean speeds
        cat(bgCyan(paste0("Speed  ")), "\n")
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = 4, trace = FALSE) #in m/s
        
        
        #weighted proportion of land use
        cat(bgGreen(paste0("Proportion of Use  ")), "\n")
        HR <- rast(raster(AKDES, DF = "PMF"))
        HR2 <- project(HR, crs(cover_2024), res = res(cover_2024))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(cover_2024), res = res(cover_2024))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        HR.df$land_class <- raster::extract(cover_2024, HR.df[,1:2])[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("9","1","3", "4","5","6","49","29")] <- "Forested"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        
        #Use the home range PDF to calculate the weighted percentage of time spent the different land class types
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
      
        
        # store results in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        #RSFs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- RSF
        use_land[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- USE
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rda for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  #save(RSFs, file = paste0(dir_path, "RSF_Fit/RSF_", DATA@info[1], ".rda"))
  save(use_land, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  
  # clean up environment
  rm(FIT, AKDES, SPEED, USE) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

# clear environment space
rm(FITS, wAKDEs, speed_mean, use_land, DATA_24)
gc()


 


