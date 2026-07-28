# 7 Extracting Moving Window Analysis Results
# This script includes the following:
#     1. Pull lists of individual results (i.e. Nayeli's models)
#     2. Pull each date range for an individual (i.e. iteration one of Nayeli)
#     3. Summarize results in an easily interpretable dataframe

# load packages
library(ctmm) # working with ctmm objects
library(dplyr) # building dataframes
library(lubridate) # dates for the moving window model
library(purrr) # binding dataframes together

# load tel data
load("./DATA/Wild_raised_tel_data/Data_telemetry.rda") #wild-raised 
load("./DATA/Orphaned_tel_data/Data_telemetry.rda") #orphaned 

# create file path for importing window analysis results
USE_path <- paste("./RESULTS/Window/Land_Use", sep = "") # use same path for distance to land type
UD_path <- paste("./RESULTS/Window/AKDEs", sep = "")
FIT_path <- paste ("./RESULTS/Window/Fits", sep = "")
SPEED_path<- paste("./RESULTS/Window/Mean_Speed", sep = "")

# wild-raised -----
# create dataframes
USE_wild <- data.frame()
SPEED_wild <- data.frame()
UD_wild <- data.frame()
FIT_wild <- data.frame()


for(i in 1:length(DATA_wild)){
  # load land use filepath
  tryCatch({
    USEPATH <- file.path(USE_path,
                            paste("Use_",
                                  DATA_wild[[i]]@info[1],
                                  ".rda",
                                  sep = ""))
    # extract an individual 
    USE <- get(load(USEPATH))
  }, error = function(e) {
    cat("Individual does not have model", i, "-", e$message, "\n")
  }) # end of trycatch
  
  # create a dataframe to hold all observations from the window analysis for one individual
  use_win <- data.frame(ID = character(length(USE))) 
  for(j in 1:length(USE)){
    # pull an individual's window
    summary_df <- USE[[j]]
    # add datetime information from window
    use_win$window_start[j] <- sub(".*_(\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2})", "\\1", names(USE)[j])
    # pull land use results from window
    if ("Agriculture" %in% colnames(summary_df)) {
      use_win[j, "Agriculture"] <- summary_df[, "Agriculture"]
    }else{
      use_win[j, "Agriculture"] <- 0
    }
    if ("Development" %in% colnames(summary_df)) {
      use_win[j, "Development"] <- summary_df[, "Development"]
    }else{
      use_win[j, "Development"] <- 0
    }
    if ("Forestry" %in% colnames(summary_df)) {
      use_win[j, "Forestry"] <- summary_df[, "Forestry"]
    }else{
      use_win[j, "Forestry"] <- 0
    }
    if ("Native_Forest" %in% colnames(summary_df)) {
      use_win[j, "Native_Forest"] <- summary_df[, "Native_Forest"]
    }else{
      use_win[j, "Native_Forest"] <- 0
    }
    if ("Pasture" %in% colnames(summary_df)) {
      use_win[j, "Pasture"] <- summary_df[, "Pasture"]
    }else{
      use_win[j, "Pasture"] <- 0
    }
    if ("Water" %in% colnames(summary_df)) {
      use_win[j, "Water"] <- summary_df[, "Water"]
    }else{
      use_win[j, "Water"] <- 0
    }
  } # end of inner loop
  use_win$ID <- DATA_wild[[i]]@info[1]
  # bind individual to total results dataframe
  USE_wild <- rbind(USE_wild, use_win) # stringAsFactors argument allows us to just keep adding the new windows to the dataframe
  
  # load AKDE path
  AKDEPATH <- file.path(UD_path,
                        paste("UD_",
                              DATA_wild[[i]]@info[1],
                              ".rda",
                              sep = ""))
  # load one individual's AKDEs
  AKDE <- get(load(AKDEPATH))
  # make a dataframe long enough to hold all windows from the window analysis
  akde_win <- data.frame(ID = character(length(AKDE))) 
  for(j in 1:length(AKDE)){
    # add a window to the dataframe
    akde_win$window_start[j] <- sub(".*_(\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2})", "\\1", names(AKDE)[j])
    # attempt to grab UD for individual
    tryCatch({
    summary_df <- data.frame(t(summary(AKDE[[j]], units = FALSE)[["CI"]])) 
    }, error = function(e) {
      cat("Individual does not have model", i, "-", e$message, "\n")
    }) # end of trycatch
    # pull information from UD
    if ("area..square.meters." %in% colnames(summary_df)) {
      akde_win[j, "HR_low"] <- summary_df[1, "area..square.meters."]
    }else{
      akde_win[j, "HR_low"] <- NA
    }
    if ("area..square.meters." %in% colnames(summary_df)) {
      akde_win[j, "HR_est"] <- summary_df[1, "area..square.meters."]
    }else{
      akde_win[j, "HR_est"] <- NA
    }
    if ("area..square.meters." %in% colnames(summary_df)) {
      akde_win[j, "HR_high"] <- summary_df[1, "area..square.meters."]
    }else{
      akde_win[j, "HR_high"] <- NA
    }
  } # end of inner loop
  # add ID to UD
  akde_win$ID <- DATA_wild[[i]]@info[1]
  # bind results together
  UD_wild<- rbind(UD_wild, akde_win)
  
  
  # load speed path
  SPEEDPATH <- file.path(SPEED_path,
                         paste("speed_",
                               DATA_wild[[i]]@info[1],
                               ".rda",
                               sep = ""))
  # load individual's speeds
  SPEED <- get(load(SPEEDPATH))
  # create a dataframe to house individual's results
  speed_win <- data.frame(ID = character(length(SPEED))) 
  for(j in 1:length(SPEED)){
    # add window information to dataframe
    speed_win$window_start[j] <- sub(".*_(\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2})", "\\1", names(SPEED)[j])
    # pull window of speed
    speedind <- SPEED[[j]]
    # extract speed results
    df <- speedind[[2]] 
    if ("low" %in% colnames(df)) {
      speed_win[j, "Mean_Speed_low"] <- df[1, "low"]
    }else{
      speed_win[j, "Mean_Speed_low"] <- NA
    }
    if ("est" %in% colnames(df)) {
      speed_win[j, "Mean_Speed_est"] <- df[1, "est"]
    }else{
      speed_win[j, "Mean_Speed_est"] <- NA
    }
    if ("high" %in% colnames(df)) {
      speed_win[j, "Mean_Speed_high"] <- df[1, "high"]
    }else{
      speed_win[j, "Mean_Speed_high"] <- NA
    }
  } # end of inner loop
  speed_win$ID <- DATA_wild[[i]]@info[1]
  # bind to total dataframe
  SPEED_wild <- rbind(SPEED_wild, speed_win)
  
  
  # load movement model path
  FITPATH <- file.path(FIT_path,
                       paste("FIT_",
                             DATA_wild[[i]]@info[1],
                             ".rda",
                             sep = ""))
  # load individual's fit
  FIT <- get(load(FITPATH))
  # make a dataframe long enough to hold all windows from the window analysis
  fit_win <- data.frame(ID = character(length(FIT))) 
  for(j in 1:length(FIT)){
    # extract window from individual
    fits <- FIT[[j]]
    # add a datetime of window to dataframe
    fit_win$window_start[j] <- sub(".*_(\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2})", "\\1", names(FIT)[j])
    # extract results of movement model for window
    summary_df <- data.frame(t(summary(fits, units = FALSE)[["CI"]])) 
    if ("area..square.meters." %in% colnames(summary_df)) {
      fit_win[j, "area_low"] <- summary_df["low", "area..square.meters."]
    }else{
      fit_win[j, "area_low"] <- NA
    }
    if ("area..square.meters." %in% colnames(summary_df)) {
      fit_win[j, "area_est"] <- summary_df["est", "area..square.meters."]
    }else{
      fit_win[j, "area_est"] <- NA
    }
    if ("area..square.meters." %in% colnames(summary_df)) {
      fit_win[j, "area_high"] <- summary_df["high", "area..square.meters."]
    }else{
      fit_win[j, "area_high"] <- NA
    }
    if ("τ.position...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τposition_low"] <- summary_df["low", "τ.position...seconds."]
    }else{
      fit_win[j, "τposition_low"] <- NA
    }
    if ("τ.position...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τposition_est"] <- summary_df["est", "τ.position...seconds."]
    }else{
      fit_win[j, "τposition_est"] <- NA
    }
    if ("τ.position...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τposition_high"] <- summary_df["high", "τ.position...seconds."]
    }else{
      fit_win[j, "τposition_high"] <- NA
    }
    if ("τ.velocity...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τvelocity_low"] <- summary_df["low", "τ.velocity...seconds."]
    }else{
      fit_win[j, "τvelocity_low"] <- NA
    }
    if ("τ.velocity...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τvelocity_est"] <- summary_df["est", "τ.velocity...seconds."]
    }else{
      fit_win[j, "τvelocity_est"] <- NA
    }
    if ("τ.velocity...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τvelocity_high"] <- summary_df["high", "τ.velocity...seconds."]
    }else{
      fit_win[j, "τvelocity_high"] <- NA
    }
    if ("speed..meters.second." %in% colnames(summary_df)) {
      fit_win[j, "speed_low"] <- summary_df["low", "speed..meters.second."]
    }else{
      fit_win[j, "speed_low"] <- NA
    } 
    if ("speed..meters.second." %in% colnames(summary_df)) {
      fit_win[j, "speed_est"] <- summary_df["est", "speed..meters.second."]
    }else{
      fit_win[j, "speed_est"] <- NA
    } 
    if ("speed..meters.second." %in% colnames(summary_df)) {
      fit_win[j, "speed_high"] <- summary_df["high", "speed..meters.second."]
    }else{
      fit_win[j, "speed_high"] <- NA
    } 
    if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
      fit_win[j, "diffusion_low"] <- summary_df["low", "diffusion..square.meters.second."]
    }else{
      fit_win[j, "diffusion_low"] <- NA
    } 
    if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
      fit_win[j, "diffusion_est"] <- summary_df["est", "diffusion..square.meters.second."]
    }else{
      fit_win[j, "diffusion_est"] <- NA
    } 
    if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
      fit_win[j, "diffusion_high"] <- summary_df["high", "diffusion..square.meters.second."]
    }else{
      fit_win[j, "diffusion_high"] <- NA
    } 
  } # end of inner loop
  fit_win$ID <- DATA_wild[[i]]@info[1]
  # bind dataframe to total results
  FIT_wild <- rbind(FIT_wild, fit_win)

  
} # end of outer loop

# create a list of all dataframes
list_df <- list(FIT_wild, UD_wild, SPEED_wild, USE_wild)
# bind the windows together
window_wild <- list_df %>%
  reduce(full_join, by = c("ID", "window_start"))

# amend columns for dates
window_wild$window_start <- as.POSIXct(window_wild$window_start, format = "%Y-%m-%d %H:%M:%S")
window_wild$month_day <- format(window_wild$window_start, "%m/%d")
window_wild$month_day <- as.Date(window_wild$month_day, format = "%m/%d")
window_wild$numeric_date <- as.numeric(window_wild$month_day) 

# convert units
# area
window_wild[,"area_low"] <- "km^2" %#% window_wild[,"area_low"] # convert m^2 to km^2 
window_wild[, "area_est"] <- "km^2" %#% window_wild[,"area_est"] # convert m^2 to km^2 
window_wild[, "area_high"] <- "km^2" %#% window_wild[,"area_high"] # convert m^2 to km^2 
# τposition
window_wild[, "τposition_low"] <- "day" %#% window_wild[,"τposition_low"] # converts seconds to days 
window_wild[, "τposition_est"] <- "day" %#% window_wild[,"τposition_est"] # converts seconds to days 
window_wild[, "τposition_high"] <- "day" %#% window_wild[,"τposition_high"] # converts seconds to days 
# τvelocity
window_wild[, "τvelocity_low"] <- "minutes" %#% window_wild[,"τvelocity_low"] # converts seconds to minutes
window_wild[, "τvelocity_est"] <- "minutes" %#% window_wild[,"τvelocity_est"] # converts seconds to minutes
window_wild[, "τvelocity_high"] <- "minutes" %#% window_wild[,"τvelocity_high"] # converts seconds to minutes
# speed
window_wild[, "speed_low"] <- "km/day" %#% window_wild[, "speed_low"] # converts m/s to km/day
window_wild[, "speed_est"] <- "km/day" %#% window_wild[, "speed_est"] # converts m/s to km/day
window_wild[, "speed_high"] <- "km/day" %#% window_wild[, "speed_high"] # converts m/s to km/day
# diffusion
window_wild[, "diffusion_low"] <- "km^2/day" %#% window_wild[, "diffusion_low"] # converts tm^2/so km^2/day 
window_wild[, "diffusion_est"] <- "km^2/day" %#% window_wild[, "diffusion_est"] # converts tm^2/so km^2/day 
window_wild[, "diffusion_high"] <- "km^2/day" %#% window_wild[, "diffusion_high"] # converts tm^2/so km^2/day 
# HR
window_wild[,"HR_low"] <- "km^2" %#% window_wild[,"HR_low"] # converts m^2 to km^2
window_wild[,"HR_est"] <- "km^2" %#% window_wild[,"HR_est"] # converts m^2 to km^2
window_wild[,"HR_high"] <- "km^2" %#% window_wild[,"HR_high"] # converts m^2 to km^2
# mean speed
window_wild[, "Mean_Speed_low"] <- "km/day" %#% window_wild[, "Mean_Speed_low"] # converts m/s to km/day
window_wild[, "Mean_Speed_est"] <- "km/day" %#% window_wild[, "Mean_Speed_est"] # converts m/s to km/day
window_wild[, "Mean_Speed_high"] <- "km/day" %#% window_wild[, "Mean_Speed_high"] # converts m/s to km/day

# save output
save(window_wild, file = "./RESULTS/Window/Window_analysis_wild_raised.rda")




# orphaned ----
USE_orphan <- data.frame()
SPEED_orphan <- data.frame()
UD_orphan <- data.frame()
FIT_orphan <- data.frame()


for(i in 1:length(DATA_orphan)){
  # load land use filepath
  tryCatch({
    USEPATH <- file.path(USE_path,
                         paste("Use_",
                               DATA_orphan[[i]]@info[1],
                               ".rda",
                               sep = ""))
    # extract an individual 
    USE <- get(load(USEPATH))
  }, error = function(e) {
    cat("Individual does not have model", i, "-", e$message, "\n")
  }) # end of trycatch
  
  # create a dataframe to hold all observations from the window analysis for one individual
  use_win <- data.frame(ID = character(length(USE))) 
  for(j in 1:length(USE)){
    # pull an individual's window
    summary_df <- USE[[j]]
    # add datetime information from window
    use_win$window_start[j] <- sub(".*_(\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2})", "\\1", names(USE)[j])
    # pull land use results from window
    if ("Agriculture" %in% colnames(summary_df)) {
      use_win[j, "Agriculture"] <- summary_df[, "Agriculture"]
    }else{
      use_win[j, "Agriculture"] <- 0
    }
    if ("Development" %in% colnames(summary_df)) {
      use_win[j, "Development"] <- summary_df[, "Development"]
    }else{
      use_win[j, "Development"] <- 0
    }
    if ("Forestry" %in% colnames(summary_df)) {
      use_win[j, "Forestry"] <- summary_df[, "Forestry"]
    }else{
      use_win[j, "Forestry"] <- 0
    }
    if ("Native_Forest" %in% colnames(summary_df)) {
      use_win[j, "Native_Forest"] <- summary_df[, "Native_Forest"]
    }else{
      use_win[j, "Native_Forest"] <- 0
    }
    if ("Pasture" %in% colnames(summary_df)) {
      use_win[j, "Pasture"] <- summary_df[, "Pasture"]
    }else{
      use_win[j, "Pasture"] <- 0
    }
    if ("Water" %in% colnames(summary_df)) {
      use_win[j, "Water"] <- summary_df[, "Water"]
    }else{
      use_win[j, "Water"] <- 0
    }
  } # end of inner loop
  use_win$ID <- DATA_orphan[[i]]@info[1]
  # bind individual to total results dataframe
  USE_orphan <- rbind(USE_orphan, use_win) # stringAsFactors argument allows us to just keep adding the new windows to the dataframe
  
  # load AKDE path
  AKDEPATH <- file.path(UD_path,
                        paste("UD_",
                              DATA_orphan[[i]]@info[1],
                              ".rda",
                              sep = ""))
  # load one individual's AKDEs
  AKDE <- get(load(AKDEPATH))
  # make a dataframe long enough to hold all windows from the window analysis
  akde_win <- data.frame(ID = character(length(AKDE))) 
  for(j in 1:length(AKDE)){
    # add a window to the dataframe
    akde_win$window_start[j] <- sub(".*_(\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2})", "\\1", names(AKDE)[j])
    # attempt to grab UD for individual
    tryCatch({
      summary_df <- data.frame(t(summary(AKDE[[j]], units = FALSE)[["CI"]])) 
    }, error = function(e) {
      cat("Individual does not have model", i, "-", e$message, "\n")
    }) # end of trycatch
    # pull information from UD
    if ("area..square.meters." %in% colnames(summary_df)) {
      akde_win[j, "HR_low"] <- summary_df[1, "area..square.meters."]
    }else{
      akde_win[j, "HR_low"] <- NA
    }
    if ("area..square.meters." %in% colnames(summary_df)) {
      akde_win[j, "HR_est"] <- summary_df[1, "area..square.meters."]
    }else{
      akde_win[j, "HR_est"] <- NA
    }
    if ("area..square.meters." %in% colnames(summary_df)) {
      akde_win[j, "HR_high"] <- summary_df[1, "area..square.meters."]
    }else{
      akde_win[j, "HR_high"] <- NA
    }
  } # end of inner loop
  # add ID to UD
  akde_win$ID <- DATA_orphan[[i]]@info[1]
  # bind results together
  UD_orphan<- rbind(UD_orphan, akde_win)
  
  
  # load speed path
  SPEEDPATH <- file.path(SPEED_path,
                         paste("speed_",
                               DATA_orphan[[i]]@info[1],
                               ".rda",
                               sep = ""))
  # load individual's speeds
  SPEED <- get(load(SPEEDPATH))
  # create a dataframe to house individual's results
  speed_win <- data.frame(ID = character(length(SPEED))) 
  for(j in 1:length(SPEED)){
    # add window information to dataframe
    speed_win$window_start[j] <- sub(".*_(\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2})", "\\1", names(SPEED)[j])
    # pull window of speed
    speedind <- SPEED[[j]]
    # extract speed results
    df <- speedind[[2]] 
    if ("low" %in% colnames(df)) {
      speed_win[j, "Mean_Speed_low"] <- df[1, "low"]
    }else{
      speed_win[j, "Mean_Speed_low"] <- NA
    }
    if ("est" %in% colnames(df)) {
      speed_win[j, "Mean_Speed_est"] <- df[1, "est"]
    }else{
      speed_win[j, "Mean_Speed_est"] <- NA
    }
    if ("high" %in% colnames(df)) {
      speed_win[j, "Mean_Speed_high"] <- df[1, "high"]
    }else{
      speed_win[j, "Mean_Speed_high"] <- NA
    }
  } # end of inner loop
  speed_win$ID <- DATA_orphan[[i]]@info[1]
  # bind to total dataframe
  SPEED_orphan <- rbind(SPEED_orphan, speed_win)
  
  
  # load movement model path
  FITPATH <- file.path(FIT_path,
                       paste("FIT_",
                             DATA_orphan[[i]]@info[1],
                             ".rda",
                             sep = ""))
  # load individual's fit
  FIT <- get(load(FITPATH))
  # make a dataframe long enough to hold all windows from the window analysis
  fit_win <- data.frame(ID = character(length(FIT))) 
  for(j in 1:length(FIT)){
    # extract window from individual
    fits <- FIT[[j]]
    # add a datetime of window to dataframe
    fit_win$window_start[j] <- sub(".*_(\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2})", "\\1", names(FIT)[j])
    # extract results of movement model for window
    summary_df <- data.frame(t(summary(fits, units = FALSE)[["CI"]])) 
    if ("area..square.meters." %in% colnames(summary_df)) {
      fit_win[j, "area_low"] <- summary_df["low", "area..square.meters."]
    }else{
      fit_win[j, "area_low"] <- NA
    }
    if ("area..square.meters." %in% colnames(summary_df)) {
      fit_win[j, "area_est"] <- summary_df["est", "area..square.meters."]
    }else{
      fit_win[j, "area_est"] <- NA
    }
    if ("area..square.meters." %in% colnames(summary_df)) {
      fit_win[j, "area_high"] <- summary_df["high", "area..square.meters."]
    }else{
      fit_win[j, "area_high"] <- NA
    }
    if ("τ.position...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τposition_low"] <- summary_df["low", "τ.position...seconds."]
    }else{
      fit_win[j, "τposition_low"] <- NA
    }
    if ("τ.position...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τposition_est"] <- summary_df["est", "τ.position...seconds."]
    }else{
      fit_win[j, "τposition_est"] <- NA
    }
    if ("τ.position...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τposition_high"] <- summary_df["high", "τ.position...seconds."]
    }else{
      fit_win[j, "τposition_high"] <- NA
    }
    if ("τ.velocity...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τvelocity_low"] <- summary_df["low", "τ.velocity...seconds."]
    }else{
      fit_win[j, "τvelocity_low"] <- NA
    }
    if ("τ.velocity...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τvelocity_est"] <- summary_df["est", "τ.velocity...seconds."]
    }else{
      fit_win[j, "τvelocity_est"] <- NA
    }
    if ("τ.velocity...seconds." %in% colnames(summary_df)) {
      fit_win[j, "τvelocity_high"] <- summary_df["high", "τ.velocity...seconds."]
    }else{
      fit_win[j, "τvelocity_high"] <- NA
    }
    if ("speed..meters.second." %in% colnames(summary_df)) {
      fit_win[j, "speed_low"] <- summary_df["low", "speed..meters.second."]
    }else{
      fit_win[j, "speed_low"] <- NA
    } 
    if ("speed..meters.second." %in% colnames(summary_df)) {
      fit_win[j, "speed_est"] <- summary_df["est", "speed..meters.second."]
    }else{
      fit_win[j, "speed_est"] <- NA
    } 
    if ("speed..meters.second." %in% colnames(summary_df)) {
      fit_win[j, "speed_high"] <- summary_df["high", "speed..meters.second."]
    }else{
      fit_win[j, "speed_high"] <- NA
    } 
    if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
      fit_win[j, "diffusion_low"] <- summary_df["low", "diffusion..square.meters.second."]
    }else{
      fit_win[j, "diffusion_low"] <- NA
    } 
    if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
      fit_win[j, "diffusion_est"] <- summary_df["est", "diffusion..square.meters.second."]
    }else{
      fit_win[j, "diffusion_est"] <- NA
    } 
    if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
      fit_win[j, "diffusion_high"] <- summary_df["high", "diffusion..square.meters.second."]
    }else{
      fit_win[j, "diffusion_high"] <- NA
    } 
  } # end of inner loop
  fit_win$ID <- DATA_orphan[[i]]@info[1]
  # bind dataframe to total results
  FIT_orphan <- rbind(FIT_orphan, fit_win)
  
  
} # end of outer loop

# create a list of all dataframes
list_df <- list(FIT_orphan, UD_orphan, SPEED_orphan, USE_orphan)
# bind the windows together
window_orphan <- list_df %>%
  reduce(full_join, by = c("ID", "window_start"))

# amend columns for dates
window_orphan$window_start <- as.POSIXct(window_orphan$window_start, format = "%Y-%m-%d %H:%M:%S")
window_orphan$month_day <- format(window_orphan$window_start, "%m/%d")
window_orphan$month_day <- as.Date(window_orphan$month_day, format = "%m/%d")
window_orphan$numeric_date <- as.numeric(window_orphan$month_day) 

# convert units
# area
window_orphan[,"area_low"] <- "km^2" %#% window_orphan[,"area_low"] # convert m^2 to km^2 
window_orphan[, "area_est"] <- "km^2" %#% window_orphan[,"area_est"] # convert m^2 to km^2 
window_orphan[, "area_high"] <- "km^2" %#% window_orphan[,"area_high"] # convert m^2 to km^2 
# τposition
window_orphan[, "τposition_low"] <- "day" %#% window_orphan[,"τposition_low"] # converts seconds to days 
window_orphan[, "τposition_est"] <- "day" %#% window_orphan[,"τposition_est"] # converts seconds to days 
window_orphan[, "τposition_high"] <- "day" %#% window_orphan[,"τposition_high"] # converts seconds to days 
# τvelocity
window_orphan[, "τvelocity_low"] <- "minutes" %#% window_orphan[,"τvelocity_low"] # converts seconds to minutes
window_orphan[, "τvelocity_est"] <- "minutes" %#% window_orphan[,"τvelocity_est"] # converts seconds to minutes
window_orphan[, "τvelocity_high"] <- "minutes" %#% window_orphan[,"τvelocity_high"] # converts seconds to minutes
# speed
window_orphan[, "speed_low"] <- "km/day" %#% window_orphan[, "speed_low"] # converts m/s to km/day
window_orphan[, "speed_est"] <- "km/day" %#% window_orphan[, "speed_est"] # converts m/s to km/day
window_orphan[, "speed_high"] <- "km/day" %#% window_orphan[, "speed_high"] # converts m/s to km/day
# diffusion
window_orphan[, "diffusion_low"] <- "km^2/day" %#% window_orphan[, "diffusion_low"] # converts tm^2/so km^2/day 
window_orphan[, "diffusion_est"] <- "km^2/day" %#% window_orphan[, "diffusion_est"] # converts tm^2/so km^2/day 
window_orphan[, "diffusion_high"] <- "km^2/day" %#% window_orphan[, "diffusion_high"] # converts tm^2/so km^2/day 
# HR
window_orphan[,"HR_low"] <- "km^2" %#% window_orphan[,"HR_low"] # converts m^2 to km^2
window_orphan[,"HR_est"] <- "km^2" %#% window_orphan[,"HR_est"] # converts m^2 to km^2
window_orphan[,"HR_high"] <- "km^2" %#% window_orphan[,"HR_high"] # converts m^2 to km^2
# mean speed
window_orphan[, "Mean_Speed_low"] <- "km/day" %#% window_orphan[, "Mean_Speed_low"] # converts m/s to km/day
window_orphan[, "Mean_Speed_est"] <- "km/day" %#% window_orphan[, "Mean_Speed_est"] # converts m/s to km/day
window_orphan[, "Mean_Speed_high"] <- "km/day" %#% window_orphan[, "Mean_Speed_high"] # converts m/s to km/day

# save output
save(window_orphan, file = "./RESULTS/Window_analysis_orphaned.rda")


