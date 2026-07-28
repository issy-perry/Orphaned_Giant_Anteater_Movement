# 8 Extract META from Moving Window Model
# This script includes the following:
#     1. Pull individuals and their windows from 
#     2. Pull one window from each individual
#     3. Calculate META from windows


# load packages
library(ctmm)
library(stringr) 

# load telemetry data for names
load("./DATA/Wild_raised_tel_data/Data_telemetry_RR.rda") # wild-raised 
load("./DATA/Orphaned_tel_data/Data_telemetry_RR.rda") # orphaned 

# load comparison data of population average
load("./RESULTS/Fits/Fits_wild_RR.rda")
load("./RESULTS/Speed/Speed_wild_RR.rda")

# remove fractal models from speed
SPEED_wild <- SPEED_wild[-c(8,26)] # Elaine and Segre

# create folder paths for pulling window analysis results
fits_folder = "./RESULTS/Window/Fits"
akde_folder = "./RESULTS/Window/AKDEs"
speed_folder = "./RESULTS/Window/Mean_Speed"

UD_path <- paste("./RESULTS/Window/AKDEs", sep = "")
FIT_path <- paste ("./RESULTS/Window/Fits", sep = "")
SPEED_path<- paste("./RESULTS/Window/Mean_Speed", sep = "")


# load windows into nested list ----
# Movement models
# get file names and their paths
fits_o <- file.path(fits_folder, paste0("FIT_", names(DATA_orphan), ".rda")) 

# pull fits for orphans
FITS_orphan <- lapply(fits_o, function(cur_file) {
  load(cur_file)
  return(FITS)
})

# rename the lists based on identity rather than another named list as these are now alphabetical
for(i in 1:length(FITS_orphan)) {
  FITS_ID <- FITS_orphan[[i]]
  names(FITS_orphan)[i] <- FITS_ID[[1]]@info$identity
}


# AKDEs
# make file path of orphans
akde_o <- file.path(akde_folder, paste0("UD_", names(DATA_orphan), ".rda")) 

# pull UDs for orphans
AKDE_orphan <- lapply(akde_o, function(cur_file) {
  load(cur_file)
  return(wAKDEs)
})

# rename the lists based on identity rather than another named list as these are now alphabetical
for(i in 1:length(AKDE_orphan)) {
  FITS_ID <- FITS_orphan[[i]]
  names(AKDE_orphan)[i] <- FITS_ID[[1]]@info$identity
}


# make file path of wild-raised individuals
akde_w <- file.path(akde_folder, paste0("UD_", names(DATA_wild), ".rda"))

# load UDs for wild-raised individuals
AKDE_wild <- lapply(akde_w, function(cur_file) {
  load(cur_file)
  return(wAKDEs)
})

# rename the lists based on identity rather than another named list as these are now alphabetical

  names(AKDE_wild) <- names(FITS_wild)


# Speeds
# make file path of orphans
speed_o <- file.path(speed_folder, paste0("Speed_", names(DATA_orphan), ".rda")) 

# pull speeds for orphans
SPEED_orphan <- lapply(speed_o, function(cur_file) {
  load(cur_file)
  return(speed_mean)
})

# rename the lists based on identity rather than another named list as these are now alphabetical
for(i in 1:length(SPEED_orphan)) {
  FITS_ID <- FITS_orphan[[i]]
  names(SPEED_orphan)[i] <- FITS_ID[[1]]@info$identity
}

# META ----
# create a dataframe to house the information
META_win_orphan <- data.frame()

# for.loop pulls windows from individuals
for(j in 1:123){ # number chosen as it is the second highest number of lists in  the nested list
  
  # make list to hold extracted windows
  fits_win <- list()
  
  # work within a new for.loop to pull windows from each individual 
  for(i in 1:length(FITS_orphan)){
    
    # pull out individual 
    IND <- FITS_orphan[[i]] 
    
    # pull out window and stick it into a slot on new list
    if(j <= length(IND)) {fits_win[[length(fits_win)+1]] <- IND[[j]]} 
  } # end inner inner loop
  
  #extract AKDEs
  akde_win <- list()
  for(i in 1:length(AKDE_orphan)){
    # pull out individual 
    IND <- AKDE_orphan[[i]] 
    
    # pull out window and stick it into a slot on new list
    if(j <= length(IND)) {akde_win[[length(akde_win)+1]] <- IND[[j]]}
  } # end inner inner loop
  
  
  akde_wild <- list()
  for(h in 1:length(AKDE_wild)){
    # pull out individual
    IND <- AKDE_wild[[h]] 
    
    # pull out window and stick it into a slot on new list
    if(j <= length(IND)) {akde_wild[[length(akde_wild)+1]] <- IND[[j]]} 
    
  } # end inner inner loop
  
  speed_win <- list()
  for(i in 1:length(SPEED_orphan)){
    # pull out individual 
    IND <- SPEED_orphan[[i]]
    
    # pull out window and stick it into a slot on new list
    if(j <= length(IND)) {speed_win[[length(speed_win)+1]] <- IND[[j]]}
  } # end inner inner loop
  
  # ensure only non fractal models are included in analysis
  akde_win <- akde_win[sapply(akde_win, function(x) length(x) >= 13)] # non fractal wAKDE object has 13 elements
  akde_wild <- akde_wild[sapply(akde_wild, function(x) length(x) >= 13)] # non fractal wAKDE object has 13 elements
  speed_win <- speed_win[sapply(speed_win, function(x) {x[[1]] != 0})]  # DOF values should not be 0
  
  
  # home range
  if (length(akde_win) > 2) {
  HR <- ctmm::meta(list(orphans = akde_win, wild = akde_wild), variable = "area", level = 0.95, units = FALSE) # units don't matter, it's a ratio
  HR_o <- as.data.frame(HR)
  HR_o$HR_low <- HR_o["orphans/", "/wild.low"]
  HR_o$HR_est <- HR_o["orphans/", "/wild.est"]
  HR_o$HR_high <- HR_o["orphans/", "/wild.high"]
  HR_o <- HR_o[-c(2), -c(1,2,3,4,5,6)] # removes column 2 and first six rows
  } else {
  HR_o <- data.frame(HR_low = NA,
                      HR_est = NA,
                      HR_high = NA)
  }
  
  # diffusion
  if (length(fits_win) > 2) {
  DIFFUSION <- ctmm::meta(list(orphans = fits_win, wild = FITS_wild), variable = "diffusion", level = 0.95, units = FALSE) # units don't matter, it's a ratio
  DIFFUSION_o <- as.data.frame(DIFFUSION)
  DIFFUSION_o$diffusion_low <- DIFFUSION_o["orphans/", "/wild.low"]
  DIFFUSION_o$diffusion_est <- DIFFUSION_o["orphans/", "/wild.est"]
  DIFFUSION_o$diffusion_high <- DIFFUSION_o["orphans/", "/wild.high"]
  DIFFUSION_o <- DIFFUSION_o[-c(2), -c(1,2,3,4,5,6)] # removes row two and first six columns 
  } else {
  DIFFUSION_o <- data.frame(diffusion_low = NA,
                            diffusion_est = NA,
                            diffusion_high= NA)
  }
  # speed
  if (length(speed_win) > 2) {
  SPEED <- ctmm::meta(list(orphans = speed_win, wild = SPEED_wild), variable = "speed", level = 0.95, units = FALSE) # units don't matter, it's a ratio
  SPEED_o <- as.data.frame(SPEED)
  SPEED_o$speed_low <- SPEED_o["orphans/", "/wild.low"]
  SPEED_o$speed_est <- SPEED_o["orphans/", "/wild.est"]
  SPEED_o$speed_high <- SPEED_o["orphans/", "/wild.high"]
  SPEED_o <- SPEED_o[-c(2), -c(1,2,3,4,5,6)] # removes row two and first six columns 
  } else {
  SPEED_o <- data.frame(speed_low = NA,
                        speed_est = NA,
                        speed_high = NA)
  }
  # tau velocity
  tryCatch({
  if (length(fits_win) > 2) {
  TAUV <- ctmm::meta(list(orphans = fits_win, wild = FITS_wild), variable = "tau velocity", level = 0.95, units = FALSE) # units don't matter, it's a ratio
  TAUV_o <- as.data.frame(TAUV)
  TAUV_o$tauvelocity_low <- TAUV_o["orphans/", "/wild.low"]
  TAUV_o$tauvelocity_est <- TAUV_o["orphans/", "/wild.est"]
  TAUV_o$tauvelocity_high <- TAUV_o["orphans/", "/wild.high"]
  TAUV_o <- TAUV_o[-c(2), -c(1,2,3,4,5,6)] # removes row two and first six columns 
  } else {
  TAUV_o <- data.frame(tauvelocity_low = NA,
                        tauvelocity_est = NA,
                        tauvelocity_high = NA)
  }
}, error = function(e) {
  cat("Error running META on Tau velocity (creating NA dataframe):", j, "-", e$message, "\n")
  TAUV_o <- data.frame(tauvelocity_low = NA,
                       tauvelocity_est = NA,
                       tauvelocity_high = NA)
})
  
  # make dataframe for window 1's results
  META_o <- cbind(HR_o, TAUV_o, DIFFUSION_o, SPEED_o)
  
  # bind to total dataframe
  META_win_orphan <- rbind(META_win_orphan, META_o) 
  
}#end of outer loop

# remove rownames
rownames(META_win_orphan) <- NULL

# remove rows with NA values for all values
META_win_orphan <- META_win_orphan[!is.na(META_win_orphan$HR_low),]

# add information regarding dNULL# add information regarding days since release
META_win_orphan$time_since_release = seq(10, by = 2, length.out = nrow(META_win_orphan)) # days since release for plotting 

# save output
save(META_win_orphan, file = "./RESULTS/Window/Total_meta_RR_df.rda")

