# 6 Window Analysis
# This script includes the following:
#     1. Extract land types from rasters to create covariates for RSFs
#     2. Fit RSFs using tel data and AKDEs against each covariate

#The moving window code is primarily adapted from Kat Chhen's analysis for her masters thesis 
#This code is based on Ryan Gill's moving window (which in turn is based on Dr. Michael Noonan's code)

#load packages 
library(ctmm) #working with ctmm telemetry objects, fitting movement models and UDs, calculating mean speeds, and using RSFs
library(lutz) #working with times
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(lubridate) #working with dates
library(tictoc) #seeing how long each individual takes to run
library(sf)
library(terra) #working with Spatrasters
library(raster) #converting Spatrasters to raster classes for RSFs
library(crayon) #adding colored bands to printed announcements for tictoc
library(weights) # For calculating the weighted proportions of the land class types


#load telemetry data
load("./DATA/Orphaned/Data_telemetry.rda") #orphaned 
load("./DATA/Wild_raised/Data_telemetry.rda") #wild-raised


#rename for convenience
DATA_wild <- DATA_TELEMETRY
rm(DATA_TELEMETRY) #clear up environment space



#individuals need to be separated by year so that 
DATA_17 <- DATA_wild[c("Anthony", "Bumpus", "Cate", "Christoffer", "Elaine", "Jackson", "Kyle", "Little_Rick", "Makao", 
                       "Puji", "Segre")]
DATA_18 <- DATA_wild[c("Alexander", "Annie", "Beto", "Hannah", "Jane", "Larry", "Luigi", "Margaret", "Maria", "Reid", 
                       "Rodolfo", "Sheron", "Thomas", "Delphine", "Gala")]
#orphaned
DATA_19 <- DATA_orphan[c("Arya", "Capitu", "Dumbo_1", "Dumbo_2")]
DATA_20 <- DATA_orphan[c("Tim_1")] 
DATA_21 <- DATA_orphan[c("Renee_1", "Renee_2", "Renee_3", "Renee_4", "Tim_2")]
DATA_22 <- DATA_orphan[c("Cláudio", "Colete", "Heather", "Juju_1", "Juju_2", "Mulan", "Peter", "Rita", "Tim_3")] 
DATA_23 <- DATA_orphan[c("Bella", "Erick", "George", "Nancy")]
DATA_24 <- DATA_orphan[c("Bahia", "Beezie", "Dom", "Jacobina", "Nayeli")] #2025 is Jacobina

#clear environment space
rm(DATA_wild, DATA_orphan)



#set time periods for window
#window was determined by META output of tau position for range-resident orphans (2.1 days est) and wild-raised(1.4 days est) 
dt <- 1 %#% "day" #shift 1 day
win <- 2 %#% "day" #window covers 2 days

#create folders to hold results
folder_list <- c("Fits", 
                 "AKDEs",
                 "Mean_Speed",
                 "Land_Use",
                 "RSF_Fit")


#set directory path
dir_path <- "./RESULTS/Windows/"


# create every folder in the folder list in the directory path
for (folder in folder_list) {
  dir.create(paste0(dir_path, folder), 
             recursive = TRUE, showWarnings = TRUE)
}



#2017  ------------------------------------------------------------------------------------------------------------------------
#set up rasters
cover_2017 <- rast("./DATA/Mapbiomas/2017_cover.tif")
CRS_rast <- cover_2017
crs(cover_2017) <- "local"
`%notin%` <- Negate(`%in%`)

Agriculture <- cover_2017
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- NA
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- 1
Agriculture <- terra::distance(Agriculture)
crs(Agriculture) <- CRS_rast
Agriculture <- raster(Agriculture)

Development <- cover_2017
Development[Development %notin% c(24,25,30)] <- NA
Development[Development %in% c(24,25,30)] <- 1
Development <- terra::distance(Development)
crs(Development) <- CRS_rast
Development <- raster(Development)

Forestry <- cover_2017
Forestry[Forestry %notin% c(9)] <- NA
Forestry[Forestry %in% c(9)] <- 1
Forestry <- terra::distance(Forestry)
crs(Forestry) <- CRS_rast
Forestry <- raster(Forestry)

Native_Forest <- cover_2017
Native_Forest[Native_Forest %notin% c(1,3,4,5,6,49,29)] <- NA
Native_Forest[Native_Forest %in% c(1,3,4,5,6,49,29)] <- 1
Native_Forest <- terra::distance(Native_Forest)
crs(Native_Forest) <- CRS_rast
Native_Forest <- raster(Native_Forest)

Pasture <- cover_2017
Pasture[Pasture %notin% c(12,15)] <- NA
Pasture[Pasture %in% c(12,15)] <- 1
Pasture <- terra::distance(Pasture)
crs(Pasture) <- CRS_rast
Pasture <- raster(Pasture)

Water <- cover_2017
Water[Water %notin% c(11,26,33)] <- NA
Water[Water %in% c(11,26,33)] <- 1
Water <- terra::distance(Water)
crs(Water) <- CRS_rast
Water <- raster(Water)

cover_2017 <- raster(cover_2017)

#window analysis
for(i in 1:length(DATA_17)){
  tic("window analysis")
  # subset out an individual
  DATA <- DATA_17[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 1 day forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), origin = "1970-01-01", tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                                                                          DATA$longitude[1],
                                                                                                          method = "accurate")) 
  
  # set up list to store
  FITS <- list()
  wAKDEs <- list()
  speed_mean <- list()
  USE <- list()
  DIST <- list()
  
  
  
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
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 3, cores = -1))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        AKDES <- akde(SUBSET, FIT, weights = TRUE)
        
        #mean speeds
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = -1) #in m/s
        
        #weighted proportion of land use
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Land_values <- raster::extract(cover_2017, HR.df[,1:2], cellnumbers = TRUE)
        HR.df$land_class <- Land_values[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        #combinind pasture and agriculture for noncover and forestry and native forest for cover
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29","9")] <- "Cover"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48","12","15")] <- "Noncover"
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        PROPS2 <- data.frame(class = names(PROPS),
                             proportion = as.numeric(PROPS))
        PROPS <- data.frame(t(PROPS2))[2,]
        names(PROPS) <- PROPS2$class
        #make dataframe
        res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        res$ID <- UD@info$identity
        #bind results to dataframe
        res <- cbind(res,PROPS)
        
        # now extract use of distance to land type
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Agriculture_values <- raster::extract(Agriculture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Agriculture <- Agriculture_values[,2]
        Development_values <- raster::extract(Development, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Development <- Development_values[,2]
        Forestry_values <- raster::extract(Forestry, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Forestry <- Forestry_values[,2]
        Native_Forest_values <- raster::extract(Native_Forest, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Native_Forest <- Native_Forest_values[,2]
        Pasture_values <- raster::extract(Pasture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Pasture <- Pasture_values[,2]
        Water_values <- raster::extract(Water, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Water <- Water_values[,2]
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        #PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        #PROPS2 <- data.frame(class = names(PROPS),
        #                     proportion = as.numeric(PROPS))
        #PROPS <- data.frame(t(PROPS2))[2,]
        #names(PROPS) <- PROPS2$class
        dist_res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        dist_res$ID <- UD@info$identity
        dist_res$Min_Agriculture <- min(HR.df2$Agriculture)
        dist_res$Mean_Agriculture <- sum(HR.df2$layer*HR.df2$Agriculture)
        dist_res$Max_Agriculture <- max(HR.df2$Agriculture)
        dist_res$Min_Development <- min(HR.df2$Development)
        dist_res$Mean_Development <- sum(HR.df2$layer*HR.df2$Development)
        dist_res$Max_Development <- max(HR.df2$Development)
        dist_res$Min_Forestry <- min(HR.df2$Forestry)
        dist_res$Mean_Forestry <- sum(HR.df2$layer*HR.df2$Forestry)
        dist_res$Max_Forestry <- max(HR.df2$Forestry)
        dist_res$Min_Native_Forest <- min(HR.df2$Native_Forest)
        dist_res$Mean_Native_Forest <- sum(HR.df2$layer*HR.df2$Native_Forest)
        dist_res$Max_Native_Forest <- max(HR.df2$Native_Forest)
        dist_res$Min_Pasture <- min(HR.df2$Pasture)
        dist_res$Mean_Pasture <- sum(HR.df2$layer*HR.df2$Pasture)
        dist_res$Max_Pasture <- max(HR.df2$Pasture)
        dist_res$Min_Water <- min(HR.df2$Water)
        dist_res$Mean_Water <- sum(HR.df2$layer*HR.df2$Water)
        dist_res$Max_Water <- max(HR.df2$Water)
        #res <- cbind(res,PROPS)
        
        
        
        
        
        # store models/UDs in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        USE[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- res
        DIST[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- dist_res
        
        
        
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rds for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  save(USE, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  save(DIST, file = paste0(dir_path, "Land_Use/Dist_", DATA@info[1], ".rda"))
  
  # clean up environment
  rm(FIT, AKDES, SPEED, RSF, dist_res, res) 

  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

rm(FITS, wAKDEs, speed_mean, USE, DIST, DATA_17, cover_2017)






#2018------------------------------------------------------------------------------------------------------------------------
#set up rasters
cover_2018 <- rast("./DATA/Mapbiomas/2018_cover.tif")
CRS_rast <- cover_2018
crs(cover_2018) <- "local"
`%notin%` <- Negate(`%in%`)

Agriculture <- cover_2018
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- NA
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- 1
Agriculture <- terra::distance(Agriculture)
crs(Agriculture) <- CRS_rast
Agriculture <- raster(Agriculture)

Development <- cover_2018
Development[Development %notin% c(24,25,30)] <- NA
Development[Development %in% c(24,25,30)] <- 1
Development <- terra::distance(Development)
crs(Development) <- CRS_rast
Development <- raster(Development)

Forestry <- cover_2018
Forestry[Forestry %notin% c(9)] <- NA
Forestry[Forestry %in% c(9)] <- 1
Forestry <- terra::distance(Forestry)
crs(Forestry) <- CRS_rast
Forestry <- raster(Forestry)

Native_Forest <- cover_2018
Native_Forest[Native_Forest %notin% c(1,3,4,5,6,49,29)] <- NA
Native_Forest[Native_Forest %in% c(1,3,4,5,6,49,29)] <- 1
Native_Forest <- terra::distance(Native_Forest)
crs(Native_Forest) <- CRS_rast
Native_Forest <- raster(Native_Forest)

Pasture <- cover_2018
Pasture[Pasture %notin% c(12,15)] <- NA
Pasture[Pasture %in% c(12,15)] <- 1
Pasture <- terra::distance(Pasture)
crs(Pasture) <- CRS_rast
Pasture <- raster(Pasture)

Water <- cover_2018
Water[Water %notin% c(11,26,33)] <- NA
Water[Water %in% c(11,26,33)] <- 1
Water <- terra::distance(Water)
crs(Water) <- CRS_rast
Water <- raster(Water)

cover_2018 <- raster(cover_2018)

#window analysis
for(i in 1:length(DATA_18)){
  tic("window analysis")
  # subset out an individual
  DATA <- DATA_18[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 1 day forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), origin = "1970-01-01", tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                                                                          DATA$longitude[1],
                                                                                                          method = "accurate")) 
  
  # set up list to store
  FITS <- list()
  wAKDEs <- list()
  speed_mean <- list()
  USE <- list()
  DIST <- list()
  
  
  
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
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 3, cores = -1))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        AKDES <- akde(SUBSET, FIT, weights = TRUE)
        
        #mean speeds
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = -1) #in m/s
        
        #weighted proportion of land use
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Land_values <- raster::extract(cover_2018, HR.df[,1:2], cellnumbers = TRUE)
        HR.df$land_class <- Land_values[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        #combinind pasture and agriculture for noncover and forestry and native forest for cover
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29","9")] <- "Cover"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48","12","15")] <- "Noncover"
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        PROPS2 <- data.frame(class = names(PROPS),
                             proportion = as.numeric(PROPS))
        PROPS <- data.frame(t(PROPS2))[2,]
        names(PROPS) <- PROPS2$class
        #make dataframe
        res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        res$ID <- UD@info$identity
        #bind results to dataframe
        res <- cbind(res,PROPS)
        
        # now extract use of distance to land type
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Agriculture_values <- raster::extract(Agriculture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Agriculture <- Agriculture_values[,2]
        Development_values <- raster::extract(Development, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Development <- Development_values[,2]
        Forestry_values <- raster::extract(Forestry, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Forestry <- Forestry_values[,2]
        Native_Forest_values <- raster::extract(Native_Forest, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Native_Forest <- Native_Forest_values[,2]
        Pasture_values <- raster::extract(Pasture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Pasture <- Pasture_values[,2]
        Water_values <- raster::extract(Water, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Water <- Water_values[,2]
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        #PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        #PROPS2 <- data.frame(class = names(PROPS),
        #                     proportion = as.numeric(PROPS))
        #PROPS <- data.frame(t(PROPS2))[2,]
        #names(PROPS) <- PROPS2$class
        dist_res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        dist_res$ID <- UD@info$identity
        dist_res$Min_Agriculture <- min(HR.df2$Agriculture)
        dist_res$Mean_Agriculture <- sum(HR.df2$layer*HR.df2$Agriculture)
        dist_res$Max_Agriculture <- max(HR.df2$Agriculture)
        dist_res$Min_Development <- min(HR.df2$Development)
        dist_res$Mean_Development <- sum(HR.df2$layer*HR.df2$Development)
        dist_res$Max_Development <- max(HR.df2$Development)
        dist_res$Min_Forestry <- min(HR.df2$Forestry)
        dist_res$Mean_Forestry <- sum(HR.df2$layer*HR.df2$Forestry)
        dist_res$Max_Forestry <- max(HR.df2$Forestry)
        dist_res$Min_Native_Forest <- min(HR.df2$Native_Forest)
        dist_res$Mean_Native_Forest <- sum(HR.df2$layer*HR.df2$Native_Forest)
        dist_res$Max_Native_Forest <- max(HR.df2$Native_Forest)
        dist_res$Min_Pasture <- min(HR.df2$Pasture)
        dist_res$Mean_Pasture <- sum(HR.df2$layer*HR.df2$Pasture)
        dist_res$Max_Pasture <- max(HR.df2$Pasture)
        dist_res$Min_Water <- min(HR.df2$Water)
        dist_res$Mean_Water <- sum(HR.df2$layer*HR.df2$Water)
        dist_res$Max_Water <- max(HR.df2$Water)
        #res <- cbind(res,PROPS)
        
        
        
        
        
        # store models/UDs in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        USE[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- res
        DIST[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- dist_res
        
        
        
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rds for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  save(USE, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  save(DIST, file = paste0(dir_path, "Land_Use/Dist_", DATA@info[1], ".rda"))
  
  
  # clean up environment
  rm(FIT, AKDES, SPEED, dist_res, res) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

rm(FITS, wAKDEs, speed_mean, USE, DIST, DATA_18, cover_2018)














#2019 ------------------------------------------------------------------------------------------------------------------------
#set up rasters
cover_2019 <- rast("./DATA/Mapbiomas/2019_cover.tif")
CRS_rast <- cover_2019
crs(cover_2019) <- "local"
`%notin%` <- Negate(`%in%`)

Agriculture <- cover_2019
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- NA
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- 1
Agriculture <- terra::distance(Agriculture)
crs(Agriculture) <- CRS_rast
Agriculture <- raster(Agriculture)

Development <- cover_2019
Development[Development %notin% c(24,25,30)] <- NA
Development[Development %in% c(24,25,30)] <- 1
Development <- terra::distance(Development)
crs(Development) <- CRS_rast
Development <- raster(Development)

Forestry <- cover_2019
Forestry[Forestry %notin% c(9)] <- NA
Forestry[Forestry %in% c(9)] <- 1
Forestry <- terra::distance(Forestry)
crs(Forestry) <- CRS_rast
Forestry <- raster(Forestry)

Native_Forest <- cover_2019
Native_Forest[Native_Forest %notin% c(1,3,4,5,6,49,29)] <- NA
Native_Forest[Native_Forest %in% c(1,3,4,5,6,49,29)] <- 1
Native_Forest <- terra::distance(Native_Forest)
crs(Native_Forest) <- CRS_rast
Native_Forest <- raster(Native_Forest)

Pasture <- cover_2019
Pasture[Pasture %notin% c(12,15)] <- NA
Pasture[Pasture %in% c(12,15)] <- 1
Pasture <- terra::distance(Pasture)
crs(Pasture) <- CRS_rast
Pasture <- raster(Pasture)

Water <- cover_2019
Water[Water %notin% c(11,26,33)] <- NA
Water[Water %in% c(11,26,33)] <- 1
Water <- terra::distance(Water)
crs(Water) <- CRS_rast
Water <- raster(Water)

cover_2019 <- raster(cover_2019)

#window analysis
for(i in 1:length(DATA_19)){
  tic("window analysis")
  # subset out an individual
  DATA <- DATA_19[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 1 day forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), origin = "1970-01-01", tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                                                                          DATA$longitude[1],
                                                                                                          method = "accurate")) 
  
  # set up list to store
  FITS <- list()
  wAKDEs <- list()
  speed_mean <- list()
  USE <- list()
  DIST <- list()
  
  
  
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
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 3, cores = -1))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        AKDES <- akde(SUBSET, FIT, weights = TRUE)
        
        #mean speeds
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = -1) #in m/s
        
        #weighted proportion of land use
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Land_values <- raster::extract(cover_2019, HR.df[,1:2], cellnumbers = TRUE)
        HR.df$land_class <- Land_values[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        #combinind pasture and agriculture for noncover and forestry and native forest for cover
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29","9")] <- "Cover"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48","12","15")] <- "Noncover"
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        PROPS2 <- data.frame(class = names(PROPS),
                             proportion = as.numeric(PROPS))
        PROPS <- data.frame(t(PROPS2))[2,]
        names(PROPS) <- PROPS2$class
        #make dataframe
        res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        res$ID <- UD@info$identity
        #bind results to dataframe
        res <- cbind(res,PROPS)
        
        # now extract use of distance to land type
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Agriculture_values <- raster::extract(Agriculture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Agriculture <- Agriculture_values[,2]
        Development_values <- raster::extract(Development, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Development <- Development_values[,2]
        Forestry_values <- raster::extract(Forestry, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Forestry <- Forestry_values[,2]
        Native_Forest_values <- raster::extract(Native_Forest, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Native_Forest <- Native_Forest_values[,2]
        Pasture_values <- raster::extract(Pasture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Pasture <- Pasture_values[,2]
        Water_values <- raster::extract(Water, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Water <- Water_values[,2]
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        #PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        #PROPS2 <- data.frame(class = names(PROPS),
        #                     proportion = as.numeric(PROPS))
        #PROPS <- data.frame(t(PROPS2))[2,]
        #names(PROPS) <- PROPS2$class
        dist_res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        dist_res$ID <- UD@info$identity
        dist_res$Min_Agriculture <- min(HR.df2$Agriculture)
        dist_res$Mean_Agriculture <- sum(HR.df2$layer*HR.df2$Agriculture)
        dist_res$Max_Agriculture <- max(HR.df2$Agriculture)
        dist_res$Min_Development <- min(HR.df2$Development)
        dist_res$Mean_Development <- sum(HR.df2$layer*HR.df2$Development)
        dist_res$Max_Development <- max(HR.df2$Development)
        dist_res$Min_Forestry <- min(HR.df2$Forestry)
        dist_res$Mean_Forestry <- sum(HR.df2$layer*HR.df2$Forestry)
        dist_res$Max_Forestry <- max(HR.df2$Forestry)
        dist_res$Min_Native_Forest <- min(HR.df2$Native_Forest)
        dist_res$Mean_Native_Forest <- sum(HR.df2$layer*HR.df2$Native_Forest)
        dist_res$Max_Native_Forest <- max(HR.df2$Native_Forest)
        dist_res$Min_Pasture <- min(HR.df2$Pasture)
        dist_res$Mean_Pasture <- sum(HR.df2$layer*HR.df2$Pasture)
        dist_res$Max_Pasture <- max(HR.df2$Pasture)
        dist_res$Min_Water <- min(HR.df2$Water)
        dist_res$Mean_Water <- sum(HR.df2$layer*HR.df2$Water)
        dist_res$Max_Water <- max(HR.df2$Water)
        #res <- cbind(res,PROPS)
        
        
        
        
        
        # store models/UDs in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        USE[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- res
        DIST[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- dist_res
        
        
        
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rds for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  save(USE, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  save(DIST, file = paste0(dir_path, "Land_Use/Dist_", DATA@info[1], ".rda"))
  
  
  # clean up environment
  rm(FIT, AKDES, SPEED, dist_res, res) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

rm(FITS, wAKDEs, speed_mean, USE, DIST, DATA_19, cover_2019)









#2020  ------------------------------------------------------------------------------------------------------------------------
#set up rasters
cover_2020 <- rast("./DATA/Mapbiomas/2020_cover.tif")
CRS_rast <- cover_2020
crs(cover_2020) <- "local"
`%notin%` <- Negate(`%in%`)

Agriculture <- cover_2020
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- NA
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- 1
Agriculture <- terra::distance(Agriculture)
crs(Agriculture) <- CRS_rast
Agriculture <- raster(Agriculture)

Development <- cover_2020
Development[Development %notin% c(24,25,30)] <- NA
Development[Development %in% c(24,25,30)] <- 1
Development <- terra::distance(Development)
crs(Development) <- CRS_rast
Development <- raster(Development)

Forestry <- cover_2020
Forestry[Forestry %notin% c(9)] <- NA
Forestry[Forestry %in% c(9)] <- 1
Forestry <- terra::distance(Forestry)
crs(Forestry) <- CRS_rast
Forestry <- raster(Forestry)

Native_Forest <- cover_2020
Native_Forest[Native_Forest %notin% c(1,3,4,5,6,49,29)] <- NA
Native_Forest[Native_Forest %in% c(1,3,4,5,6,49,29)] <- 1
Native_Forest <- terra::distance(Native_Forest)
crs(Native_Forest) <- CRS_rast
Native_Forest <- raster(Native_Forest)

Pasture <- cover_2020
Pasture[Pasture %notin% c(12,15)] <- NA
Pasture[Pasture %in% c(12,15)] <- 1
Pasture <- terra::distance(Pasture)
crs(Pasture) <- CRS_rast
Pasture <- raster(Pasture)

Water <- cover_2020
Water[Water %notin% c(11,26,33)] <- NA
Water[Water %in% c(11,26,33)] <- 1
Water <- terra::distance(Water)
crs(Water) <- CRS_rast
Water <- raster(Water)

cover_2020 <- raster(cover_2020)

#window analysis
for(i in 1:length(DATA_20)){
  tic("window analysis")
  # subset out an individual
  DATA <- DATA_20[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 1 day forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), origin = "1970-01-01", tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                                                                          DATA$longitude[1],
                                                                                                          method = "accurate")) 
  
  # set up list to store
  FITS <- list()
  wAKDEs <- list()
  speed_mean <- list()
  USE <- list()
  DIST <- list()
  
  
  
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
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 3, cores = -1))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        AKDES <- akde(SUBSET, FIT, weights = TRUE)
        
        #mean speeds
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = -1) #in m/s
        
        #weighted proportion of land use
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Land_values <- raster::extract(cover_2020, HR.df[,1:2], cellnumbers = TRUE)
        HR.df$land_class <- Land_values[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        #combinind pasture and agriculture for noncover and forestry and native forest for cover
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29","9")] <- "Cover"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48","12","15")] <- "Noncover"
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        PROPS2 <- data.frame(class = names(PROPS),
                             proportion = as.numeric(PROPS))
        PROPS <- data.frame(t(PROPS2))[2,]
        names(PROPS) <- PROPS2$class
        #make dataframe
        res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        res$ID <- UD@info$identity
        #bind results to dataframe
        res <- cbind(res,PROPS)
        
        # now extract use of distance to land type
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Agriculture_values <- raster::extract(Agriculture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Agriculture <- Agriculture_values[,2]
        Development_values <- raster::extract(Development, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Development <- Development_values[,2]
        Forestry_values <- raster::extract(Forestry, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Forestry <- Forestry_values[,2]
        Native_Forest_values <- raster::extract(Native_Forest, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Native_Forest <- Native_Forest_values[,2]
        Pasture_values <- raster::extract(Pasture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Pasture <- Pasture_values[,2]
        Water_values <- raster::extract(Water, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Water <- Water_values[,2]
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        #PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        #PROPS2 <- data.frame(class = names(PROPS),
        #                     proportion = as.numeric(PROPS))
        #PROPS <- data.frame(t(PROPS2))[2,]
        #names(PROPS) <- PROPS2$class
        dist_res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        dist_res$ID <- UD@info$identity
        dist_res$Min_Agriculture <- min(HR.df2$Agriculture)
        dist_res$Mean_Agriculture <- sum(HR.df2$layer*HR.df2$Agriculture)
        dist_res$Max_Agriculture <- max(HR.df2$Agriculture)
        dist_res$Min_Development <- min(HR.df2$Development)
        dist_res$Mean_Development <- sum(HR.df2$layer*HR.df2$Development)
        dist_res$Max_Development <- max(HR.df2$Development)
        dist_res$Min_Forestry <- min(HR.df2$Forestry)
        dist_res$Mean_Forestry <- sum(HR.df2$layer*HR.df2$Forestry)
        dist_res$Max_Forestry <- max(HR.df2$Forestry)
        dist_res$Min_Native_Forest <- min(HR.df2$Native_Forest)
        dist_res$Mean_Native_Forest <- sum(HR.df2$layer*HR.df2$Native_Forest)
        dist_res$Max_Native_Forest <- max(HR.df2$Native_Forest)
        dist_res$Min_Pasture <- min(HR.df2$Pasture)
        dist_res$Mean_Pasture <- sum(HR.df2$layer*HR.df2$Pasture)
        dist_res$Max_Pasture <- max(HR.df2$Pasture)
        dist_res$Min_Water <- min(HR.df2$Water)
        dist_res$Mean_Water <- sum(HR.df2$layer*HR.df2$Water)
        dist_res$Max_Water <- max(HR.df2$Water)
        #res <- cbind(res,PROPS)
        
        
        
        
        
        # store models/UDs in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        USE[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- res
        DIST[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- dist_res
        
        
        
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rds for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  save(USE, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  save(DIST, file = paste0(dir_path, "Land_Use/Dist_", DATA@info[1], ".rda"))
  
  # clean up environment
  rm(FIT, AKDES, SPEED, dist_res, res) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

rm(FITS, wAKDEs, speed_mean, USE, DIST, DATA_20, cover_2020)








#2021  ------------------------------------------------------------------------------------------------------------------------
#set up rasters
cover_2021 <- rast("./DATA/Mapbiomas/2021_cover.tif")
CRS_rast <- cover_2021
crs(cover_2021) <- "local"
`%notin%` <- Negate(`%in%`)

Agriculture <- cover_2021
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- NA
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- 1
Agriculture <- terra::distance(Agriculture)
crs(Agriculture) <- CRS_rast
Agriculture <- raster(Agriculture)

Development <- cover_2021
Development[Development %notin% c(24,25,30)] <- NA
Development[Development %in% c(24,25,30)] <- 1
Development <- terra::distance(Development)
crs(Development) <- CRS_rast
Development <- raster(Development)

Forestry <- cover_2021
Forestry[Forestry %notin% c(9)] <- NA
Forestry[Forestry %in% c(9)] <- 1
Forestry <- terra::distance(Forestry)
crs(Forestry) <- CRS_rast
Forestry <- raster(Forestry)

Native_Forest <- cover_2021
Native_Forest[Native_Forest %notin% c(1,3,4,5,6,49,29)] <- NA
Native_Forest[Native_Forest %in% c(1,3,4,5,6,49,29)] <- 1
Native_Forest <- terra::distance(Native_Forest)
crs(Native_Forest) <- CRS_rast
Native_Forest <- raster(Native_Forest)

Pasture <- cover_2021
Pasture[Pasture %notin% c(12,15)] <- NA
Pasture[Pasture %in% c(12,15)] <- 1
Pasture <- terra::distance(Pasture)
crs(Pasture) <- CRS_rast
Pasture <- raster(Pasture)

Water <- cover_2021
Water[Water %notin% c(11,26,33)] <- NA
Water[Water %in% c(11,26,33)] <- 1
Water <- terra::distance(Water)
crs(Water) <- CRS_rast
Water <- raster(Water)

cover_2021 <- raster(cover_2021)

#window analysis
for(i in 1:length(DATA_21)){
  tic("window analysis")
  # subset out an individual
  DATA <- DATA_21[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 1 day forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), origin = "1970-01-01", tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                                                                          DATA$longitude[1],
                                                                                                          method = "accurate")) 
  
  # set up list to store
  FITS <- list()
  wAKDEs <- list()
  speed_mean <- list()
  USE <- list()
  DIST <- list()
  
  
  
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
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 3, cores = -1))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        AKDES <- akde(SUBSET, FIT, weights = TRUE)
        
        #mean speeds
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = -1) #in m/s
        
        #weighted proportion of land use
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Land_values <- raster::extract(cover_2021, HR.df[,1:2], cellnumbers = TRUE)
        HR.df$land_class <- Land_values[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        #combinind pasture and agriculture for noncover and forestry and native forest for cover
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29","9")] <- "Cover"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48","12","15")] <- "Noncover"
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        PROPS2 <- data.frame(class = names(PROPS),
                             proportion = as.numeric(PROPS))
        PROPS <- data.frame(t(PROPS2))[2,]
        names(PROPS) <- PROPS2$class
        #make dataframe
        res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        res$ID <- UD@info$identity
        #bind results to dataframe
        res <- cbind(res,PROPS)
        
        # now extract use of distance to land type
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Agriculture_values <- raster::extract(Agriculture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Agriculture <- Agriculture_values[,2]
        Development_values <- raster::extract(Development, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Development <- Development_values[,2]
        Forestry_values <- raster::extract(Forestry, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Forestry <- Forestry_values[,2]
        Native_Forest_values <- raster::extract(Native_Forest, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Native_Forest <- Native_Forest_values[,2]
        Pasture_values <- raster::extract(Pasture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Pasture <- Pasture_values[,2]
        Water_values <- raster::extract(Water, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Water <- Water_values[,2]
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        #PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        #PROPS2 <- data.frame(class = names(PROPS),
        #                     proportion = as.numeric(PROPS))
        #PROPS <- data.frame(t(PROPS2))[2,]
        #names(PROPS) <- PROPS2$class
        dist_res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        dist_res$ID <- UD@info$identity
        dist_res$Min_Agriculture <- min(HR.df2$Agriculture)
        dist_res$Mean_Agriculture <- sum(HR.df2$layer*HR.df2$Agriculture)
        dist_res$Max_Agriculture <- max(HR.df2$Agriculture)
        dist_res$Min_Development <- min(HR.df2$Development)
        dist_res$Mean_Development <- sum(HR.df2$layer*HR.df2$Development)
        dist_res$Max_Development <- max(HR.df2$Development)
        dist_res$Min_Forestry <- min(HR.df2$Forestry)
        dist_res$Mean_Forestry <- sum(HR.df2$layer*HR.df2$Forestry)
        dist_res$Max_Forestry <- max(HR.df2$Forestry)
        dist_res$Min_Native_Forest <- min(HR.df2$Native_Forest)
        dist_res$Mean_Native_Forest <- sum(HR.df2$layer*HR.df2$Native_Forest)
        dist_res$Max_Native_Forest <- max(HR.df2$Native_Forest)
        dist_res$Min_Pasture <- min(HR.df2$Pasture)
        dist_res$Mean_Pasture <- sum(HR.df2$layer*HR.df2$Pasture)
        dist_res$Max_Pasture <- max(HR.df2$Pasture)
        dist_res$Min_Water <- min(HR.df2$Water)
        dist_res$Mean_Water <- sum(HR.df2$layer*HR.df2$Water)
        dist_res$Max_Water <- max(HR.df2$Water)
        #res <- cbind(res,PROPS)
        
        
        
        
        
        # store models/UDs in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        USE[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- res
        DIST[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- dist_res
        
        
        
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rds for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  save(USE, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  save(DIST, file = paste0(dir_path, "Land_Use/Dist_", DATA@info[1], ".rda"))
  
  
  # clean up environment
  rm(FIT, AKDES, SPEED, dist_res, res) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

rm(FITS, wAKDEs, speed_mean, USE, DIST, DATA_21, cover_2021)










#2022  ------------------------------------------------------------------------------------------------------------------------
#set up rasters
cover_2022 <- rast("./DATA/Mapbiomas/2022_cover.tif")
CRS_rast <- cover_2022
crs(cover_2022) <- "local"
`%notin%` <- Negate(`%in%`)

Agriculture <- cover_2022
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- NA
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- 1
Agriculture <- terra::distance(Agriculture)
crs(Agriculture) <- CRS_rast
Agriculture <- raster(Agriculture)

Development <- cover_2022
Development[Development %notin% c(24,25,30)] <- NA
Development[Development %in% c(24,25,30)] <- 1
Development <- terra::distance(Development)
crs(Development) <- CRS_rast
Development <- raster(Development)

Forestry <- cover_2022
Forestry[Forestry %notin% c(9)] <- NA
Forestry[Forestry %in% c(9)] <- 1
Forestry <- terra::distance(Forestry)
crs(Forestry) <- CRS_rast
Forestry <- raster(Forestry)

Native_Forest <- cover_2022
Native_Forest[Native_Forest %notin% c(1,3,4,5,6,49,29)] <- NA
Native_Forest[Native_Forest %in% c(1,3,4,5,6,49,29)] <- 1
Native_Forest <- terra::distance(Native_Forest)
crs(Native_Forest) <- CRS_rast
Native_Forest <- raster(Native_Forest)

Pasture <- cover_2022
Pasture[Pasture %notin% c(12,15)] <- NA
Pasture[Pasture %in% c(12,15)] <- 1
Pasture <- terra::distance(Pasture)
crs(Pasture) <- CRS_rast
Pasture <- raster(Pasture)

Water <- cover_2022
Water[Water %notin% c(11,26,33)] <- NA
Water[Water %in% c(11,26,33)] <- 1
Water <- terra::distance(Water)
crs(Water) <- CRS_rast
Water <- raster(Water)

cover_2022 <- raster(cover_2022)

#window analysis
for(i in 1:length(DATA_22)){
  tic("window analysis")
  # subset out an individual
  DATA <- DATA_22[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 1 day forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), origin = "1970-01-01", tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                                                                          DATA$longitude[1],
                                                                                                          method = "accurate")) 
  
  # set up list to store
  FITS <- list()
  wAKDEs <- list()
  speed_mean <- list()
  USE <- list()
  DIST <- list()
  
  
  
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
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 3, cores = -1))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        AKDES <- akde(SUBSET, FIT, weights = TRUE)
        
        #mean speeds
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = -1) #in m/s
        
        #weighted proportion of land use
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Land_values <- raster::extract(cover_2022, HR.df[,1:2], cellnumbers = TRUE)
        HR.df$land_class <- Land_values[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        #combinind pasture and agriculture for noncover and forestry and native forest for cover
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29","9")] <- "Cover"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48","12","15")] <- "Noncover"
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        PROPS2 <- data.frame(class = names(PROPS),
                             proportion = as.numeric(PROPS))
        PROPS <- data.frame(t(PROPS2))[2,]
        names(PROPS) <- PROPS2$class
        #make dataframe
        res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        res$ID <- UD@info$identity
        #bind results to dataframe
        res <- cbind(res,PROPS)
        
        # now extract use of distance to land type
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Agriculture_values <- raster::extract(Agriculture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Agriculture <- Agriculture_values[,2]
        Development_values <- raster::extract(Development, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Development <- Development_values[,2]
        Forestry_values <- raster::extract(Forestry, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Forestry <- Forestry_values[,2]
        Native_Forest_values <- raster::extract(Native_Forest, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Native_Forest <- Native_Forest_values[,2]
        Pasture_values <- raster::extract(Pasture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Pasture <- Pasture_values[,2]
        Water_values <- raster::extract(Water, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Water <- Water_values[,2]
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        #PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        #PROPS2 <- data.frame(class = names(PROPS),
        #                     proportion = as.numeric(PROPS))
        #PROPS <- data.frame(t(PROPS2))[2,]
        #names(PROPS) <- PROPS2$class
        dist_res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        dist_res$ID <- UD@info$identity
        dist_res$Min_Agriculture <- min(HR.df2$Agriculture)
        dist_res$Mean_Agriculture <- sum(HR.df2$layer*HR.df2$Agriculture)
        dist_res$Max_Agriculture <- max(HR.df2$Agriculture)
        dist_res$Min_Development <- min(HR.df2$Development)
        dist_res$Mean_Development <- sum(HR.df2$layer*HR.df2$Development)
        dist_res$Max_Development <- max(HR.df2$Development)
        dist_res$Min_Forestry <- min(HR.df2$Forestry)
        dist_res$Mean_Forestry <- sum(HR.df2$layer*HR.df2$Forestry)
        dist_res$Max_Forestry <- max(HR.df2$Forestry)
        dist_res$Min_Native_Forest <- min(HR.df2$Native_Forest)
        dist_res$Mean_Native_Forest <- sum(HR.df2$layer*HR.df2$Native_Forest)
        dist_res$Max_Native_Forest <- max(HR.df2$Native_Forest)
        dist_res$Min_Pasture <- min(HR.df2$Pasture)
        dist_res$Mean_Pasture <- sum(HR.df2$layer*HR.df2$Pasture)
        dist_res$Max_Pasture <- max(HR.df2$Pasture)
        dist_res$Min_Water <- min(HR.df2$Water)
        dist_res$Mean_Water <- sum(HR.df2$layer*HR.df2$Water)
        dist_res$Max_Water <- max(HR.df2$Water)
        #res <- cbind(res,PROPS)
        
        
        
        
        
        # store models/UDs in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        USE[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- res
        DIST[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- dist_res
        
        
        
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rds for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  save(USE, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  save(DIST, file = paste0(dir_path, "Land_Use/Dist_", DATA@info[1], ".rda"))
  
  
  # clean up environment
  rm(FIT, AKDES, SPEED, dist_res, res) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

rm(FITS, wAKDEs, speed_mean, USE, DIST, DATA_22, cover_2022)







#2023  ------------------------------------------------------------------------------------------------------------------------
#set up rasters
cover_2023 <- rast("./DATA/Mapbiomas/2023_cover.tif")
CRS_rast <- cover_2023
crs(cover_2023) <- "local"
`%notin%` <- Negate(`%in%`)

Agriculture <- cover_2023
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- NA
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- 1
Agriculture <- terra::distance(Agriculture)
crs(Agriculture) <- CRS_rast
Agriculture <- raster(Agriculture)

Development <- cover_2023
Development[Development %notin% c(24,25,30)] <- NA
Development[Development %in% c(24,25,30)] <- 1
Development <- terra::distance(Development)
crs(Development) <- CRS_rast
Development <- raster(Development)

Forestry <- cover_2023
Forestry[Forestry %notin% c(9)] <- NA
Forestry[Forestry %in% c(9)] <- 1
Forestry <- terra::distance(Forestry)
crs(Forestry) <- CRS_rast
Forestry <- raster(Forestry)

Native_Forest <- cover_2023
Native_Forest[Native_Forest %notin% c(1,3,4,5,6,49,29)] <- NA
Native_Forest[Native_Forest %in% c(1,3,4,5,6,49,29)] <- 1
Native_Forest <- terra::distance(Native_Forest)
crs(Native_Forest) <- CRS_rast
Native_Forest <- raster(Native_Forest)

Pasture <- cover_2023
Pasture[Pasture %notin% c(12,15)] <- NA
Pasture[Pasture %in% c(12,15)] <- 1
Pasture <- terra::distance(Pasture)
crs(Pasture) <- CRS_rast
Pasture <- raster(Pasture)

Water <- cover_2023
Water[Water %notin% c(11,26,33)] <- NA
Water[Water %in% c(11,26,33)] <- 1
Water <- terra::distance(Water)
crs(Water) <- CRS_rast
Water <- raster(Water)

cover_2023 <- raster(cover_2023)

#window analysis
for(i in 1:length(DATA_23)){
  tic("window analysis")
  # subset out an individual
  DATA <- DATA_23[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 1 day forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), origin = "1970-01-01", tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                                                                          DATA$longitude[1],
                                                                                                          method = "accurate")) 
  
  # set up list to store
  FITS <- list()
  wAKDEs <- list()
  speed_mean <- list()
  USE <- list()
  DIST <- list()
  
  
  
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
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 3, cores = -1))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        AKDES <- akde(SUBSET, FIT, weights = TRUE)
        
        #mean speeds
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = -1) #in m/s
        
        #weighted proportion of land use
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Land_values <- raster::extract(cover_2023, HR.df[,1:2], cellnumbers = TRUE)
        HR.df$land_class <- Land_values[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        #combinind pasture and agriculture for noncover and forestry and native forest for cover
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29","9")] <- "Cover"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48","12","15")] <- "Noncover"
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        PROPS2 <- data.frame(class = names(PROPS),
                             proportion = as.numeric(PROPS))
        PROPS <- data.frame(t(PROPS2))[2,]
        names(PROPS) <- PROPS2$class
        #make dataframe
        res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        res$ID <- UD@info$identity
        #bind results to dataframe
        res <- cbind(res,PROPS)
        
        # now extract use of distance to land type
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Agriculture_values <- raster::extract(Agriculture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Agriculture <- Agriculture_values[,2]
        Development_values <- raster::extract(Development, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Development <- Development_values[,2]
        Forestry_values <- raster::extract(Forestry, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Forestry <- Forestry_values[,2]
        Native_Forest_values <- raster::extract(Native_Forest, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Native_Forest <- Native_Forest_values[,2]
        Pasture_values <- raster::extract(Pasture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Pasture <- Pasture_values[,2]
        Water_values <- raster::extract(Water, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Water <- Water_values[,2]
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        #PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        #PROPS2 <- data.frame(class = names(PROPS),
        #                     proportion = as.numeric(PROPS))
        #PROPS <- data.frame(t(PROPS2))[2,]
        #names(PROPS) <- PROPS2$class
        dist_res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        dist_res$ID <- UD@info$identity
        dist_res$Min_Agriculture <- min(HR.df2$Agriculture)
        dist_res$Mean_Agriculture <- sum(HR.df2$layer*HR.df2$Agriculture)
        dist_res$Max_Agriculture <- max(HR.df2$Agriculture)
        dist_res$Min_Development <- min(HR.df2$Development)
        dist_res$Mean_Development <- sum(HR.df2$layer*HR.df2$Development)
        dist_res$Max_Development <- max(HR.df2$Development)
        dist_res$Min_Forestry <- min(HR.df2$Forestry)
        dist_res$Mean_Forestry <- sum(HR.df2$layer*HR.df2$Forestry)
        dist_res$Max_Forestry <- max(HR.df2$Forestry)
        dist_res$Min_Native_Forest <- min(HR.df2$Native_Forest)
        dist_res$Mean_Native_Forest <- sum(HR.df2$layer*HR.df2$Native_Forest)
        dist_res$Max_Native_Forest <- max(HR.df2$Native_Forest)
        dist_res$Min_Pasture <- min(HR.df2$Pasture)
        dist_res$Mean_Pasture <- sum(HR.df2$layer*HR.df2$Pasture)
        dist_res$Max_Pasture <- max(HR.df2$Pasture)
        dist_res$Min_Water <- min(HR.df2$Water)
        dist_res$Mean_Water <- sum(HR.df2$layer*HR.df2$Water)
        dist_res$Max_Water <- max(HR.df2$Water)
        #res <- cbind(res,PROPS)
        
        
        
        
        
        # store models/UDs in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        USE[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- res
        DIST[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- dist_res
        
        
        
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rds for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  save(USE, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  save(DIST, file = paste0(dir_path, "Land_Use/Dist_", DATA@info[1], ".rda"))
  
  # clean up environment
  rm(FIT, AKDES, SPEED, dist_res, res) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

rm(FITS, wAKDEs, speed_mean, USE, DIST, DATA_23, cover_2023)











#2024  ------------------------------------------------------------------------------------------------------------------------
#set up rasters
cover_2024 <- rast("./DATA/Mapbiomas/2024_cover.tif")
CRS_rast <- cover_2024
crs(cover_2024) <- "local"
`%notin%` <- Negate(`%in%`)

Agriculture <- cover_2024
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- NA
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- 1
Agriculture <- terra::distance(Agriculture)
crs(Agriculture) <- CRS_rast
Agriculture <- raster(Agriculture)

Development <- cover_2024
Development[Development %notin% c(24,25,30)] <- NA
Development[Development %in% c(24,25,30)] <- 1
Development <- terra::distance(Development)
crs(Development) <- CRS_rast
Development <- raster(Development)

Forestry <- cover_2024
Forestry[Forestry %notin% c(9)] <- NA
Forestry[Forestry %in% c(9)] <- 1
Forestry <- terra::distance(Forestry)
crs(Forestry) <- CRS_rast
Forestry <- raster(Forestry)

Native_Forest <- cover_2024
Native_Forest[Native_Forest %notin% c(1,3,4,5,6,49,29)] <- NA
Native_Forest[Native_Forest %in% c(1,3,4,5,6,49,29)] <- 1
Native_Forest <- terra::distance(Native_Forest)
crs(Native_Forest) <- CRS_rast
Native_Forest <- raster(Native_Forest)

Pasture <- cover_2024
Pasture[Pasture %notin% c(12,15)] <- NA
Pasture[Pasture %in% c(12,15)] <- 1
Pasture <- terra::distance(Pasture)
crs(Pasture) <- CRS_rast
Pasture <- raster(Pasture)

Water <- cover_2024
Water[Water %notin% c(11,26,33)] <- NA
Water[Water %in% c(11,26,33)] <- 1
Water <- terra::distance(Water)
crs(Water) <- CRS_rast
Water <- raster(Water)

cover_2024 <- raster(cover_2024)

#window analysis
for(i in 1:length(DATA_24)){
  tic("window analysis")
  # subset out an individual
  DATA <- DATA_24[[i]]
  
  # set up the window segments ......................................................
  # generate start times with a 2 day segment for the individual
  times <- seq(from = DATA$t[1], # t = Unix timestamp format
               to = DATA$t[nrow(DATA)],  
               by = dt) # shift each segment by 1 day forward
  
  #ensure that they are full days from 00:00 to 23:59, set the timestamps to 00:00 time, since we are looking at 2 day windows and not time specific
  #convert Unix timestamps to POSIXct
  times <- as.POSIXct(as.numeric(as.character(times)), origin = "1970-01-01", tz = lutz::tz_lookup_coords(DATA$latitude[1], 
                                                                                                          DATA$longitude[1],
                                                                                                          method = "accurate")) 
  
  # set up list to store
  FITS <- list()
  wAKDEs <- list()
  speed_mean <- list()
  USE <- list()
  DIST <- list()
  
  
  
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
      GUESS <- ctmm.guess(SUBSET, CTMM=ctmm(error = TRUE), interactive = FALSE)
      #movement model fits
      FIT <- try(ctmm.select(SUBSET, GUESS, trace = 3, cores = -1))
      
      if (inherits(FIT, "ctmm")) {
        #UDs
        AKDES <- akde(SUBSET, FIT, weights = TRUE)
        
        #mean speeds
        SPEED <- speed(object = SUBSET, CTMM = FIT, robust = TRUE, units = FALSE, cores = -1) #in m/s
        
        #weighted proportion of land use
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Land_values <- raster::extract(cover_2024, HR.df[,1:2], cellnumbers = TRUE)
        HR.df$land_class <- Land_values[,2]
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
        HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
        HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
        HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
        HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
        #combinind pasture and agriculture for noncover and forestry and native forest for cover
        HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29","9")] <- "Cover"
        HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48","12","15")] <- "Noncover"
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        PROPS2 <- data.frame(class = names(PROPS),
                             proportion = as.numeric(PROPS))
        PROPS <- data.frame(t(PROPS2))[2,]
        names(PROPS) <- PROPS2$class
        #make dataframe
        res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        res$ID <- UD@info$identity
        #bind results to dataframe
        res <- cbind(res,PROPS)
        
        # now extract use of distance to land type
        HR <- rast(raster(UD, DF = "PMF"))
        HR2 <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df2 <- terra::as.data.frame(HR2, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df2$layer <- HR.df2$layer/sum(HR.df2$layer)
        HR <- project(HR, crs(CRS_rast), res = res(CRS_rast))
        HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
        #Renormalize
        HR.df$layer <- HR.df$layer/sum(HR.df$layer)
        #Extract habitat values
        Agriculture_values <- raster::extract(Agriculture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Agriculture <- Agriculture_values[,2]
        Development_values <- raster::extract(Development, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Development <- Development_values[,2]
        Forestry_values <- raster::extract(Forestry, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Forestry <- Forestry_values[,2]
        Native_Forest_values <- raster::extract(Native_Forest, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Native_Forest <- Native_Forest_values[,2]
        Pasture_values <- raster::extract(Pasture, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Pasture <- Pasture_values[,2]
        Water_values <- raster::extract(Water, HR.df2[,1:2], cellnumbers = TRUE)/1000 # delete "/1000" to rescale HFI
        HR.df2$Water <- Water_values[,2]
        # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
        #PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
        #PROPS2 <- data.frame(class = names(PROPS),
        #                     proportion = as.numeric(PROPS))
        #PROPS <- data.frame(t(PROPS2))[2,]
        #names(PROPS) <- PROPS2$class
        dist_res <- data.frame(binomial = "Myrmecophaga_tridactyla")
        dist_res$ID <- UD@info$identity
        dist_res$Min_Agriculture <- min(HR.df2$Agriculture)
        dist_res$Mean_Agriculture <- sum(HR.df2$layer*HR.df2$Agriculture)
        dist_res$Max_Agriculture <- max(HR.df2$Agriculture)
        dist_res$Min_Development <- min(HR.df2$Development)
        dist_res$Mean_Development <- sum(HR.df2$layer*HR.df2$Development)
        dist_res$Max_Development <- max(HR.df2$Development)
        dist_res$Min_Forestry <- min(HR.df2$Forestry)
        dist_res$Mean_Forestry <- sum(HR.df2$layer*HR.df2$Forestry)
        dist_res$Max_Forestry <- max(HR.df2$Forestry)
        dist_res$Min_Native_Forest <- min(HR.df2$Native_Forest)
        dist_res$Mean_Native_Forest <- sum(HR.df2$layer*HR.df2$Native_Forest)
        dist_res$Max_Native_Forest <- max(HR.df2$Native_Forest)
        dist_res$Min_Pasture <- min(HR.df2$Pasture)
        dist_res$Mean_Pasture <- sum(HR.df2$layer*HR.df2$Pasture)
        dist_res$Max_Pasture <- max(HR.df2$Pasture)
        dist_res$Min_Water <- min(HR.df2$Water)
        dist_res$Mean_Water <- sum(HR.df2$layer*HR.df2$Water)
        dist_res$Max_Water <- max(HR.df2$Water)
        #res <- cbind(res,PROPS)
        
        
        
        
        
        # store models/UDs in a list, name the entry based on anteater name and subset window start date, not the times[i] as that is in unix format
        FITS[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- FIT
        wAKDEs[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- AKDES
        speed_mean[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- SPEED
        USE[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- res
        DIST[[paste0(DATA@info[1], "_", as.character(WINDOW_START))]] <- dist_res
        
        
        
        
      }
    }, error = function(e) {
      cat("Error during processing for window segment:", j, "-", e$message, "\n")
    })
    
    
  } #end of inner loop
  
  # save all the outputs as a rds for future analysis 
  message(underline(bgGreen(white(paste("saving output for anteater", DATA@info[1])))))
  save(FITS, file = paste0(dir_path, "Fits/Fit_", DATA@info[1], ".rda"))
  save(wAKDEs, file = paste0(dir_path, "AKDEs/UD_", DATA@info[1], ".rda"))
  save(speed_mean, file = paste0(dir_path, "Mean_Speed/speed_", DATA@info[1], ".rda"))
  save(USE, file = paste0(dir_path, "Land_Use/Use_", DATA@info[1], ".rda"))
  save(DIST, file = paste0(dir_path, "Land_Use/Dist_", DATA@info[1], ".rda"))
  
  
  # clean up environment
  rm(FIT, AKDES, SPEED, dist_res, res) 
  
  gc() # free up computational resources
  
  # end of outer loop, goes back to a new anteater
  toc()
  
}

rm(FITS, wAKDEs, speed_mean, USE, DIST, DATA_24, cover_2024)










































 


