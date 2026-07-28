# 1 Fit Movement Models
# This script includes the following:
#     1. Fit continuous-time movement models (hereafter ctmm) (orphaned and wild-raised populations separately)
#     2. Estimate utilization distributions, or UDs (orphaned and wild-raised populations separately)
#     3. Estimate mean speeds (orphaned and wild-raised populations separately)


# load packages
library(ctmm) # working with tel data and fitting models
library(tictoc) # determine how long each model takes to run


# load data
load("./DATA/Orphaned_tel_data/Data_telemetry.rda") # orphaned
load("./DATA/Wild_raised_tel_data/Data_telemetry.rda") # wild-raised

# orphaned----
# make lists to house results
FITS_orphan <- list()
AKDE_orphan <- list()
SPEED_orphan <- list()

# for.loop for running models for each individual
for(i in 1:length(DATA_orphan)){
  
  tic("individual")
  
  # extract individual
  DATA <- DATA_orphan[[i]]
  
  # create variograms based on individuals
  GUESS <- ctmm.guess(DATA, CTMM=ctmm(error = TRUE), interactive = FALSE) 
  
  # fit models to variograms
  FIT <- ctmm.select(DATA, GUESS, trace = TRUE, cores = -1)
  
  # add to list
  FITS_orphan[[i]] <- FIT
  
  # calculate wAKDE
  AKDE_orphan[[i]] <- akde(DATA, FIT, weights = TRUE) # weights = TRUE because we do have dispersing populations
  
  # calculate average speed and input into list
  SPEED_orphan[[i]] <- speed(object = DATA, CTMM = FIT, robust = TRUE, units = FALSE, cores = -1) #in m/s
  
  toc() # end of timing
}

#transfer names from tel data list to the movement model list
names(FITS_orphan) <- names(DATA_orphan)
names(AKDE_orphan) <- names(DATA_orphan)
names(SPEED_orphan) <- names(FITS_orphan)

#save outputs
save(FITS_orphan, file = "./RESULTS/Fits/Fits_orphan.rda") 
save(AKDE_orphan, file = "./RESULTS/AKDEs/UDs_orphan.rda")
save(SPEED_orphan, file = "./RESULTS/Speed/Speed_orphans.rda")



# wild-raised ----
# ctmm and UDs were already generated for the wild-raised population using the same methodology as above, so we just need to calculate speed here
load("./RESULTS/Fits/Fits_wild.rda") # wild-raised

# for.loop for calculating mean speed
SPEED_wild <- list()
for(i in 1:length(DATA_wild)){
  
  # extract individual's telemetry data
  DATA <- DATA_wild[[i]]
  
  # extract individual's movement model fit
  FIT <- FITS_wild[[i]]
  
  # calculate average speed and input into list
  SPEED_wild[[i]] <- speed(object = DATA, CTMM = FIT, robust = TRUE) 
}

# transfer names from tel data
names(SPEED_wild) <- names(FITS_wild)

# save output
save(SPEED_wild, file = "./RESULTS/Speed/Speed_wild_raised.rda")
