# 4 Resource Selection Functions
# This script includes the following:
#     1. Extract land types from rasters to create covariates for RSFs
#     2. Fit RSFs using tel data and AKDEs against each covariate


# load packages 
library(ctmm) # includes functions for movement models, UDs, and RSFs
library(terra) # extracting land types from Spatrasters prior to converting to raster classes
library(raster) # convert Spatrasters to raster classes for RSFs
library(lutz) #working with times
library(sf)
library(crayon) #adding colored bands to printed announcements for tictoc

# load data
load("./DATA/Wild_raised_tel_data/Data_telemetry_RR.rda") #wild-raised 
load("./DATA/Orphaned_tel_data/Data_telemetry_RR.rda") #orphaned 

#load UDs
load("./RESULTS/AKDEs/UDs_orphan_RR.rda") #orphaned 
load("./RESULTS/AKDEs/UDs_wild_RR.rda") #wild-raised 


#need to separate the data depending on when the individual was monitored so that the results for land selection are as accurate as possible
#wild-raised
DATA_17 <- DATA_wild[c("Anthony", "Bumpus", "Cate", "Christoffer", "Elaine", "Jackson", "Kyle", "Little_Rick", "Makao", "Puji")] # rm Segre disperser 
DATA_18 <- DATA_wild[c("Alexander", "Annie", "Beto", "Hannah", "Jane", "Larry", "Luigi", "Margaret", "Maria", "Reid", "Rodolfo", "Sheron", "Thomas")] # rm Gala, Delphine dispersers "Alexander", "Annie", "Beto", "Hannah", "Jane", "Larry", "Luigi", "Margaret", "Maria",
#orphaned
DATA_19 <- DATA_orphan[c("Arya", "Dumbo_1", "Dumbo_2")] # rm Capitu
#DATA_20 <- DATA_orphan[c("Tim_1")] 
DATA_21 <- DATA_orphan[c("Tim_2")] # rm "Renee_1", "Renee_2", "Renee_3", "Renee_4",
DATA_22 <- DATA_orphan[c("Colete", "Heather", "Juju_2", "Mulan", "Peter", "Rita", "Tim_3")] # rm "Cláudio","Juju_1",
DATA_23 <- DATA_orphan[c("Bella", "George", "Nancy")] # rm "Erick","George",
DATA_24 <- DATA_orphan[c("Dom")] #rm "Bahia","Beezie", "Jacobina", "Nayeli"

#subset AKDEs
#wild-raised
AKDE_17 <- AKDE_wild[c("Anthony", "Bumpus", "Cate", "Christoffer", "Elaine", "Jackson", "Kyle", "Little_Rick", "Makao", "Puji")] # rm Segre disperser 
AKDE_18 <- AKDE_wild[c( "Reid","Rodolfo", "Sheron", "Thomas")] # rm Gala, Delphine dispersers"Alexander", "Annie", "Beto", "Hannah", "Jane", "Larry", "Luigi", "Margaret", "Maria",
#orphaned
AKDE_19 <- AKDE_orphan[c("Arya", "Dumbo_1", "Dumbo_2")] # rm Capitu
#AKDE_20 <- AKDE_orphan[c("Tim_1")] 
AKDE_21 <- AKDE_orphan[c("Tim_2")] # rm "Renee_1", "Renee_2", "Renee_3", "Renee_4",
AKDE_22 <- AKDE_orphan[c("Colete", "Heather", "Juju_2", "Mulan", "Peter", "Rita", "Tim_3")] # rm "Cláudio","Juju_1",
AKDE_23 <-AKDE_orphan[c("Bella", "George", "Nancy")] # rm "Erick","George",
AKDE_24 <- AKDE_orphan[c("Dom")] #rm "Bahia","Beezie", "Jacobina", "Nayeli"
# clear up environment space
rm(DATA_wild, DATA_orphan, AKDE_orphan, AKDE_wild) 

# ensure R provides a warning before the loop is done (so that we know what individuals are problematic (if any are) so time is not wasted)
options(warn = 1)

# 2017 ---- 
# create a list to house results
RSF_17 <- list()
# for.loop for running RSFs
for(i in 1:length(DATA_17)){
  #extract individual
  DATA <- DATA_17[[i]]
  
  #extract AKDE
  AKDE <- AKDE_17[[i]]
  
  #ensure projections are the same
  ctmm::projection(DATA) <- ctmm::projection(AKDE)
  
  # keep track of individual
  cat(bgMagenta(paste("RSF for", DATA@info[1]), "\n"))
  
  #fit RSFs
  IND <- rsf.select(DATA, UD = AKDE, R = covers_2017,
                    error = 0.055, cores = 1, max.mem = "1.5 Gb")
  RSF_17[[i]] <- IND
  
  # save individual results in case of crash
  save(IND, file = paste0("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2017_", DATA@info[1], ".rda"))
  
} # close the loop

# transfer names to new RSF list
names(RSF_17) <- names(DATA_17)

# save output
save(RSF_17, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2017.rda") 

# free environment space
rm(RSF_17, DATA_17, AKDE_17, covers, cover_2017)
gc()


# 2018 ----
# create a list to house results
RSF_18 <- list()
# for.loop for running RSFs
for(i in 1:length(DATA_18)){
  #extract individual
  DATA <- DATA_18[[i]]
  
  #extract AKDE
  AKDE <- AKDE_18[[i]]
  
  #ensure projections are the same
  ctmm::projection(DATA) <- ctmm::projection(AKDE)
  
  # keep track of individual
  cat(bgMagenta(paste("RSF for", DATA@info[1]), "\n"))
  
  #fit RSF
  IND <- rsf.select(DATA, UD = AKDE, R = covers_2018, 
                    error = 0.055, cores = 1, max.mem = "1.5 Gb")
  RSF_18[[i]] <- IND
  
  # save individual results in case of crash
  save(IND, file = paste0("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2018_", DATA@info[1], ".rda"))
  
} #close the loop

#transfer names to new RSF list
names(RSF_18) <- names(DATA_18)

#save output
save(RSF_18, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2018.rda")

#free up environment space
rm(RSF_18, DATA_18, AKDE_18, covers, cover_2018, CRS_rast) 
gc()


# 2019 ----
# create a list to house results
RSF_19 <- list()
# for.loop for running RSFs
for(i in 1:length(DATA_19)){
  #extract individual
  DATA <- DATA_19[[i]]
  
  #extract AKDE
  AKDE <- AKDE_19[[i]]
  
  #ensure projections are the same
  ctmm::projection(DATA) <- ctmm::projection(AKDE)
  
  # keep track of individual
  cat(bgMagenta(paste("RSF for", DATA@info[1]), "\n"))
  
  #fit RSF
  IND <- rsf.select(DATA, UD = AKDE, R = covers_2019, 
                    error = 0.055, cores = 1, max.mem = "1.5 Gb")
  RSF_19[[i]] <- IND
  
  # save individual results in case of crash
  save(IND, file = paste0("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2019_", DATA@info[1], ".rda"))
  
} #close the loop

#transfer names to new RSF list
names(RSF_19) <- names(DATA_19)

#save output
save(RSF_19, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2019.rda") 
#free environment space
rm(RSF_19, DATA_19, AKDE_19, covers, cover_2019, CRS_rast) 
gc()


# 2021 ----
# create a list to house results
RSF_21 <- list()
# for.loop for running RSFs
for(i in 1:length(DATA_21)){
  #extract individual
  DATA <- DATA_21[[i]]
  
  #extract AKDE
  AKDE <- AKDE_21[[i]]
  
  #ensure projections are the same
  ctmm::projection(DATA) <- ctmm::projection(AKDE)
  
  # keep track of individual
  cat(bgMagenta(paste("RSF for", DATA@info[1]), "\n"))
  
  #fit RSF
  IND <- rsf.select(DATA, UD = AKDE, R = covers_2021, 
                    error = 0.055, cores = 1, max.mem = "1.5 Gb")
  RSF_21[[i]] <- IND
  
  # save individual results in case of crash
  save(IND, file = paste0("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2021_", DATA@info[1], ".rda"))
  
} #close the loop

#transfer names to new RSF list
names(RSF_21) <- names(DATA_21)

#save output
save(RSF_21, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2021.rda")
#free environment space
rm(RSF_21, DATA_21, AKDE_21, covers, cover_2021, CRS_rast) 
gc()


# 2022 ----
# create a list to house results
RSF_22 <- list()
# for.loop for running RSFs
for(i in 1:length(DATA_22)){
  #extract individual
  DATA <- DATA_22[[i]]
  
  #extract AKDE
  AKDE <- AKDE_22[[i]]
  
  #ensure projections are the same
  ctmm::projection(DATA) <- ctmm::projection(AKDE)
  
  # keep track of individual
  cat(bgMagenta(paste("RSF for", DATA@info[1]), "\n"))
  
  #fit RSF
  IND <- rsf.select(DATA, UD = AKDE, R = covers_2022, 
                    error = 0.055, cores = 1, max.mem = "1.5 Gb")
  RSF_22[[i]] <- IND
  
  # save individual results in case of crash
  save(IND, file = paste0("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2022_", DATA@info[1], ".rda"))
  
} #close the loop

#transfer names to new RSF list
names(RSF_22) <- names(DATA_22)

#save output
save(RSF_22, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2022.rda") 

#free environment space
rm(RSF_23, DATA_23, AKDE_23, covers, CRS_rast, cover_2023)
gc()


# 2023 ----
# create a list to house results
RSF_23 <- list()
# for.loop for running RSFs
for(i in 1:length(DATA_23)){
  #extract individual
  DATA <- DATA_23[[i]]
  
  #extract AKDE
  AKDE <- AKDE_23[[i]]
  
  #ensure projections are the same
  ctmm::projection(DATA) <- ctmm::projection(AKDE)
  
  # keep track of individual
  cat(bgMagenta(paste("RSF for", DATA@info[1]), "\n"))
  
  #fit RSF
  IND <- rsf.select(DATA, UD = AKDE, R = covers_2023, 
                    error = 0.055, cores = 1, max.mem = "1.5 Gb")
  RSF_23[[i]] <- IND
  
  # save individual results in case of crash
  save(IND, file = paste0("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2023_", DATA@info[1], ".rda"))

} #close the loop

#transfer names to new RSF list
names(RSF_23) <- names(DATA_23)

#save output
save(RSF_23, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2023.rda") 

#free environment space
rm(RSF_23, DATA_23, AKDE_23, covers, CRS_rast, cover_2023)
gc()


# 2024 ----
# create a list to house results
RSF_24 <- list()
# for.loop for running RSFs
for(i in 1:length(DATA_24)){
  #extract individual
  DATA <- DATA_24[[i]]
  
  #extract AKDE
  AKDE <- AKDE_24[[i]]
  
  #ensure projections are the same
  ctmm::projection(DATA) <- ctmm::projection(AKDE)
  
  # keep track of individual
  cat(bgMagenta(paste("RSF for", DATA@info[1]), "\n"))
  
  #fit RSF
  IND <- rsf.select(DATA, UD = AKDE, R = covers_2024, 
                    error = 0.055, cores = 1, max.mem = "1.5 Gb")
  RSF_24[[i]] <- IND
  
  # save individual results in case of crash
  save(IND, file = paste0("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2024_", DATA@info[1], ".rda"))
  
} #close the loop

#transfer names to new RSF list
names(RSF_24) <- names(DATA_24)
#save output
save(RSF_24, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2024.rda") 
#free environment space
rm(DATA_24, AKDE_24, covers, CRS_rast, cover_2024)
gc()


# set warnings to default again
options(warn = 0)


