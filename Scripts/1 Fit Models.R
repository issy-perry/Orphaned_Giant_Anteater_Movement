# 1 Fit Movement Models
# This script includes the following:
#     1. Fit movement models (orphaned and wild-raised populations separately)
#     2. Estimate utilization distributions, or UDs (orphaned and wild-raised populations separately)
#     3. Estimate mean speeds (orphaned and wild-raised populations separately)


# load packages
library(ctmm) # working with tel data and fitting models
library(tictoc) # prints out time a model takes for each individual


# load data
load("./DATA/Orphaned_tel_data/Data_telemetry.rda") # orphaned
load("./DATA/Wild_raised_tel_data/Data_telemetry.rda") # wild-raised
# rename for convenience 
DATA_wild <- DATA_TELEMETRY
rm(DATA_TELEMETRY)



# orphaned----
# make a list to hold the movement models
FITS_orphan <- list()
# fit movement model to each orphan
for(i in 1:length(DATA_orphan)){
  
  tic("individual")
  
  # extract individual
  DATA <- DATA_orphan[[i]]
  
  # create variograms based on individuals
  GUESS <- ctmm.guess(DATA, CTMM=ctmm(error = TRUE), interactive = FALSE) 
  
  # fit models to variograms
  FITS_orphan[[i]] <- ctmm.select(DATA, GUESS, trace = TRUE, cores = -1)
  
  toc()
}

#transfer names from tel data list to the movement model list
names(FITS_orphan) <- names(DATA_orphan)

#save output
save(FITS_orphan, file = "./RESULTS/Fits/Fits_orphan.rda") 


# Utilization distributions
# typically AKDEs are all fit in one line; however, these models are being finicky (likely due to dispersing behaviors). 
# So, we will fit individuals separately in a for.loop
# make a list to hold the UDs
AKDE_orphan <- list()

# for.loop for calculating UDs
for(i in 1:length(DATA_orphan)){
  
  # extract individual
  DATA <- DATA_orphan[[i]]
  
  # extract movement model fits
  FIT <- FITS_orphan[[i]]
  
  # calculate wAKDE
  AKDE_orphan[[i]] <- akde(DATA, FIT, weights = TRUE) # weights = TRUE because we do have dispersing populations
}

# transfer names from tel daa list to models
names(AKDE_orphan) <- names(DATA_orphan)

# save output
save(AKDE_orphan, file = "./RESULTS/AKDEs/UDs_orphan.rda")



# mean speed
# make a list to hold results
SPEED_orphan <- list()

# for.loop for calculating mean speed
for(i in 1:length(DATA_orphan)){
  
  # extract individual's telemetry data
  DATA <- DATA_orphan[[i]]
  
  # extract individual's movement model fit
  FIT <- FITS_orphan[[i]]
  
  # calculate average speed and input into list
  SPEED_orphan[[i]] <- speed(object = DATA, CTMM = FIT, robust = TRUE, units = FALSE, cores = -1) #in m/s
}

# transfer names from tel data to speed estimations
names(SPEED_orphan) <- names(FITS_orphan)

# save output
save(SPEED_orphan, file = "./RESULTS/Speed/Speed_orphans.rda")









# wild-raised ----
# movement models and UDs were already generated for the wild-raised population using the same methodology as above, so we just need to calculate speed here
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

