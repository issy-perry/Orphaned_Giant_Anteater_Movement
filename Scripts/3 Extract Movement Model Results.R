# 3 Extract Movement Model Results
# This script includes the following:
#     1. Extract results from each movement model, UD, and speed estimation
#     2. Change units so all results are comparable and easy to read
#     3. Conduct a meta-analysis on the means for each population

# load packages
library(ctmm) # working with ctmm objects

# load movement models
load("./RESULTS/Fits/Fits_orphan.rda") # orphaned 
load("./RESULTS/Fits/Fits_wild.rda") # wild-raised

# load UDs
load("./RESULTS/AKDEs/UDs_orphan.rda") # orphaned 
load("./RESULTS/AKDEs/UDs_wild.rda") # wild-raised

# load speed estimations
load("./RESULTS/Speed/Speed_orphans.rda") # orphaned 
load("./RESULTS/Speed/Speed_wild_raised.rda") # wild-raised





# orphaned population ----
# create dataframes to hold results for whole population 
FITS_orphan_df <- data.frame()
AKDE_orphan_df <- data.frame()
SPEED_orphan_df <- data.frame()

# for.loop for extracting results
for(i in 1:length(FITS_orphan)){
  # movement models ----
   # create a dataframe to hold summary outputs of movement models for EACH individual
  fit_test <- data.frame()
  
  # extract information from summary output of an individual's movement model
  summary_df <- data.frame(t(summary(FITS_orphan[[i]], units = FALSE)[["CI"]])) # units = FALSE makes all values the same scientific standard
  
  # if value is available in summary output, extract it into the individual's dataframe. 
  # if value is not available in summary output, include an NA in the column
  # area used
  if ("area..square.meters." %in% colnames(summary_df)) {
    fit_test[i, "area_low"] <- summary_df["low", "area..square.meters."]
  }else{
    fit_test[i, "area_low"] <- NA
  }
  if ("area..square.meters." %in% colnames(summary_df)) {
    fit_test[i, "area_est"] <- summary_df["est", "area..square.meters."]
  }else{
    fit_test[i, "area_est"] <- NA
  }
  if ("area..square.meters." %in% colnames(summary_df)) {
    fit_test[i, "area_high"] <- summary_df["high", "area..square.meters."]
  }else{
    fit_test[i, "area_high"] <- NA
  }
  # tau position
  if ("τ.position...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τposition_low"] <- summary_df["low", "τ.position...seconds."]
  }else{
    fit_test[i, "τposition_low"] <- NA
  }
  if ("τ.position...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τposition_est"] <- summary_df["est", "τ.position...seconds."]
  }else{
    fit_test[i, "τposition_est"] <- NA
  }
  if ("τ.position...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τposition_high"] <- summary_df["high", "τ.position...seconds."]
  }else{
    fit_test[i, "τposition_high"] <- NA
  }
  # tau velocity
  if ("τ.velocity...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τvelocity_low"] <- summary_df["low", "τ.velocity...seconds."]
  }else{
    fit_test[i, "τvelocity_low"] <- NA
  }
  if ("τ.velocity...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τvelocity_est"] <- summary_df["est", "τ.velocity...seconds."]
  }else{
    fit_test[i, "τvelocity_est"] <- NA
  }
  if ("τ.velocity...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τvelocity_high"] <- summary_df["high", "τ.velocity...seconds."]
  }else{
    fit_test[i, "τvelocity_high"] <- NA
  }
  # diffusion
  if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
    fit_test[i, "diffusion_low"] <- summary_df["low", "diffusion..square.meters.second."]
  }else{
    fit_test[i, "diffusion_low"] <- NA
  }
  if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
    fit_test[i, "diffusion_est"] <- summary_df["est", "diffusion..square.meters.second."]
  }else{
    fit_test[i, "diffusion_est"] <- NA
  }
  if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
    fit_test[i, "diffusion_high"] <- summary_df["high", "diffusion..square.meters.second."]
  }else{
    fit_test[i, "diffusion_high"] <- NA
  }
   # end of if/else loops
  
  
  # bind individuals dataframe to total dataframe
  FITS_orphan_df <- rbind(FITS_orphan_df, fit_test) 
  
  
  
  
  # UDs ----
  # create a dataframe to hold summary outputs of AKDEs for EACH individual
  akde_test <- data.frame()
  
  # extract information from summary output of an individual's AKDE
  summary_df <- data.frame(t(summary(AKDE_orphan[[i]], units = FALSE)[["CI"]])) 
  
  # if value is available in summary output, extract it into the individual's dataframe. 
  # if value is not available in summary output, include an NA in the column
  if ("area..square.meters." %in% colnames(summary_df)) {
    akde_test[i, "HR_low"] <- summary_df["low", "area..square.meters."]
  }else{
    akde_test[i, "HR_low"] <- NA
  }
  if ("area..square.meters." %in% colnames(summary_df)) {
    akde_test[i, "HR_est"] <- summary_df["est", "area..square.meters."]
  }else{
    akde_test[i, "HR_est"] <- NA
  }
  if ("area..square.meters." %in% colnames(summary_df)) {
    akde_test[i, "HR_high"] <- summary_df["high", "area..square.meters."]
  }else{
    akde_test[i, "HR_high"] <- NA
  }
  # end of if/else loops
  
  # assign name of giant anteater to the individual's dataframe
  akde_test$ID <- AKDE_orphan[[i]]@info[1]
  
  # bind individual's dataframe to total dataframe
  AKDE_orphan_df <- rbind(AKDE_orphan_df, akde_test) 
  
  
  
  
  # speed ----
  # create a dataframe to hold outputs of mean speed for EACH individual
  speed_test <- data.frame()
  
  # extract individual
  IND <- SPEED_orphan[[i]]
  
  # extract confidence intervals of speed from individual
  df <- IND[[2]] 
  
  # if value is available in summary output, extract it into the individual's dataframe. 
  # if value is not available in summary output, include an NA in the column
  if ("low" %in% colnames(df)) {
    speed_test[i, "Mean_Speed_low"] <- df[1, "low"]
  }else{
    speed_test[i, "Mean_Speed_low"] <- NA
  }
  if ("est" %in% colnames(df)) {
    speed_test[i, "Mean_Speed_est"] <- df[1, "est"]
  }else{
    speed_test[i, "Mean_Speed_est"] <- NA
  }
  if ("high" %in% colnames(df)) {
    speed_test[i, "Mean_Speed_high"] <- df[1, "high"]
  }else{
    speed_test[i, "Mean_Speed_high"] <- NA
  }
  # end of if/else loop
  
  # bind individual's dataframe to total dataframe
  SPEED_orphan_df <- rbind(SPEED_orphan_df, speed_test)

} # end of outer loop


#combine dataframes into one
Orphan_df <- cbind(FITS_orphan_df, SPEED_orphan_df, AKDE_orphan_df)

# remove extra NA values
# we will use the area column from movement models, for no matter what kind of movement model it is, it will always include area
Orphan_df <- Orphan_df[!is.na(Orphan_df$area_low),]


# convert units
# area
Orphan_df[,"area_low"] <- "km^2" %#% Orphan_df[,"area_low"] # converts m^2 to km^2 
Orphan_df[, "area_est"] <- "km^2" %#% Orphan_df[,"area_est"] # converts m^2 to km^2 
Orphan_df[, "area_high"] <- "km^2" %#% Orphan_df[,"area_high"] # converts m^2 to km^2 
# τposition
Orphan_df[, "τposition_low"] <- "day" %#% Orphan_df[,"τposition_low"] # converts seconds to days 
Orphan_df[, "τposition_est"] <- "day" %#% Orphan_df[,"τposition_est"] # converts seconds to days 
Orphan_df[, "τposition_high"] <- "day" %#% Orphan_df[,"τposition_high"] # converts seconds to days 
# τvelocity
Orphan_df[, "τvelocity_low"] <- "minutes" %#% Orphan_df[,"τvelocity_low"] # converts seconds to minutes
Orphan_df[, "τvelocity_est"] <- "minutes" %#% Orphan_df[,"τvelocity_est"] # converts seconds to minutes
Orphan_df[, "τvelocity_high"] <- "minutes" %#% Orphan_df[,"τvelocity_high"] # converts seconds to minutes
# diffusion
Orphan_df[, "diffusion_low"] <- "km^2/day" %#% Orphan_df[, "diffusion_low"] # converts m^2/s to km^2/day 
Orphan_df[, "diffusion_est"] <- "km^2/day" %#% Orphan_df[, "diffusion_est"] # converts m^2/s to km^2/day 
Orphan_df[, "diffusion_high"] <- "km^2/day" %#% Orphan_df[, "diffusion_high"] # converts m^2/s to km^2/day 
# home range size estimation
Orphan_df[,"HR_low"] <- "km^2" %#% Orphan_df[,"HR_low"] # converts m^2 to km^2 
Orphan_df[, "HR_est"] <- "km^2" %#% Orphan_df[,"HR_est"] # converts m^2 to km^2 
Orphan_df[, "HR_high"] <- "km^2" %#% Orphan_df[,"HR_high"] # converts m^2 to km^2 
# mean speed
Orphan_df[, "Mean_Speed_low"] <- "km/day" %#% Orphan_df[, "Mean_Speed_low"] # converts m/s to km/day 
Orphan_df[, "Mean_Speed_est"] <- "km/day" %#% Orphan_df[, "Mean_Speed_est"] # converts m/s to km/day 
Orphan_df[, "Mean_Speed_high"] <- "km/day" %#% Orphan_df[, "Mean_Speed_high"] # converts m/s to km/day 

# remove odd numbers for rownames
row.names(Orphan_df) <- NULL











# wild-raised population ----
# create a dataframe to hold results for whole population 
FITS_wild_df <- data.frame()
AKDE_wild_df <- data.frame()
SPEED_wild_df <- data.frame()

for(i in 1:length(FITS_wild)){
  # create a dataframe to hold summary outputs of movement models for EACH individual
  fit_test <- data.frame()
  # extract information from summary output of an individual's movement model
  summary_df <- data.frame(t(summary(FITS_wild[[i]], units = FALSE)[["CI"]]))
  # if value is available in summary output, extract it into the individual's dataframe. 
  # if value is not available in summary output, include an NA in the column
  # area used
  if ("area..square.meters." %in% colnames(summary_df)) {
    fit_test[i, "area_low"] <- summary_df["low", "area..square.meters."]
  }else{
    fit_test[i, "area_low"] <- NA
  }
  if ("area..square.meters." %in% colnames(summary_df)) {
    fit_test[i, "area_est"] <- summary_df["est", "area..square.meters."]
  }else{
    fit_test[i, "area_est"] <- NA
  }
  if ("area..square.meters." %in% colnames(summary_df)) {
    fit_test[i, "area_high"] <- summary_df["high", "area..square.meters."]
  }else{
    fit_test[i, "area_high"] <- NA
  }
  # tau position
  if ("τ.position...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τposition_low"] <- summary_df["low", "τ.position...seconds."]
  }else{
    fit_test[i, "τposition_low"] <- NA
  }
  if ("τ.position...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τposition_est"] <- summary_df["est", "τ.position...seconds."]
  }else{
    fit_test[i, "τposition_est"] <- NA
  }
  if ("τ.position...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τposition_high"] <- summary_df["high", "τ.position...seconds."]
  }else{
    fit_test[i, "τposition_high"] <- NA
  }
  # tau velocity
  if ("τ.velocity...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τvelocity_low"] <- summary_df["low", "τ.velocity...seconds."]
  }else{
    fit_test[i, "τvelocity_low"] <- NA
  }
  if ("τ.velocity...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τvelocity_est"] <- summary_df["est", "τ.velocity...seconds."]
  }else{
    fit_test[i, "τvelocity_est"] <- NA
  }
  if ("τ.velocity...seconds." %in% colnames(summary_df)) {
    fit_test[i, "τvelocity_high"] <- summary_df["high", "τ.velocity...seconds."]
  }else{
    fit_test[i, "τvelocity_high"] <- NA
  }
  #diffusion
  if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
    fit_test[i, "diffusion_low"] <- summary_df["low", "diffusion..square.meters.second."]
  }else{
    fit_test[i, "diffusion_low"] <- NA
  }
  if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
    fit_test[i, "diffusion_est"] <- summary_df["est", "diffusion..square.meters.second."]
  }else{
    fit_test[i, "diffusion_est"] <- NA
  }
  if ("diffusion..square.meters.second." %in% colnames(summary_df)) {
    fit_test[i, "diffusion_high"] <- summary_df["high", "diffusion..square.meters.second."]
  }else{
    fit_test[i, "diffusion_high"] <- NA
  }
  # end of if/else loops
  
  # bind individuals dataframe to total dataframe
  FITS_wild_df <- rbind(FITS_wild_df, fit_test) 
  
  
  # UDs ----
  # create a dataframe to hold summary outputs of AKDEs for EACH individual
  akde_test <- data.frame()
  
  # extract information from summary output of an individual's AKDE
  summary_df <- data.frame(t(summary(AKDE_wild[[i]], units = FALSE)[["CI"]])) 
  # if value is available in summary output, extract it into the individual's dataframe. 
  # if value is not available in summary output, include an NA in the column
  if ("area..square.meters." %in% colnames(summary_df)) {
    akde_test[i, "HR_low"] <- summary_df["low", "area..square.meters."]
  }else{
    akde_test[i, "HR_low"] <- NA
  }
  if ("area..square.meters." %in% colnames(summary_df)) {
    akde_test[i, "HR_est"] <- summary_df["est", "area..square.meters."]
  }else{
    akde_test[i, "HR_est"] <- NA
  }
  if ("area..square.meters." %in% colnames(summary_df)) {
    akde_test[i, "HR_high"] <- summary_df["high", "area..square.meters."]
  }else{
    akde_test[i, "HR_high"] <- NA
  }
  # end of if/else loops
  
  # assign name of giant anteater to the individual's dataframe
  akde_test$ID <- AKDE_wild[[i]]@info[1]
  
  # bind individual's dataframe to total dataframe
  AKDE_wild_df <- rbind(AKDE_wild_df, akde_test) 
  
  
  
  
  
  # mean speed ----
  # create a dataframe to house the individual's results
  speed_test <- data.frame()
  
  # extract individual
  IND <- SPEED_wild[[i]]
  
  # extract confidence intervals of speed from individual
  df <- IND[[2]] 
  # if value is available in summary output, extract it into the individual's dataframe. 
  # if value is not available in summary output, include an NA in the column
  if ("low" %in% colnames(df)) {
    speed_test[i, "Mean_Speed_low"] <- df[1, "low"]
  }else{
    speed_test[i, "Mean_Speed_low"] <- NA
  }
  if ("est" %in% colnames(df)) {
    speed_test[i, "Mean_Speed_est"] <- df[1, "est"]
  }else{
    speed_test[i, "Mean_Speed_est"] <- NA
  }
  if ("high" %in% colnames(df)) {
    speed_test[i, "Mean_Speed_high"] <- df[1, "high"]
  }else{
    speed_test[i, "Mean_Speed_high"] <- NA
  }
  # end of if/else loop
  
  # bind individual's dataframe to total dataframe
  SPEED_wild_df <- rbind(SPEED_wild_df, speed_test)
  
} # end of outer loop


# combine dataframes into one
Wild_df <- cbind(FITS_wild_df, SPEED_wild_df, AKDE_wild_df)

# remove extra NA values
Wild_df <- Wild_df[!is.na(Wild_df$area_low),]

# convert units 
# area
Wild_df[,"area_low"] <- "km^2" %#% Wild_df[,"area_low"] # converts m^2 to km^2 
Wild_df[, "area_est"] <- "km^2" %#% Wild_df[,"area_est"] # converts m^2 to km^2 
Wild_df[, "area_high"] <- "km^2" %#% Wild_df[,"area_high"] # converts m^2 to km^2 
# τposition
Wild_df[, "τposition_low"] <- "day" %#% Wild_df[,"τposition_low"] # converts seconds to days
Wild_df[, "τposition_est"] <- "day" %#% Wild_df[,"τposition_est"] # converts seconds to days
Wild_df[, "τposition_high"] <- "day" %#% Wild_df[,"τposition_high"] # converts seconds to days
# τvelocity
Wild_df[, "τvelocity_low"] <- "minutes" %#% Wild_df[,"τvelocity_low"] # converts seconds to minutes
Wild_df[, "τvelocity_est"] <- "minutes" %#% Wild_df[,"τvelocity_est"] # converts seconds to minutes
Wild_df[, "τvelocity_high"] <- "minutes" %#% Wild_df[,"τvelocity_high"] # converts seconds to minutes
# diffusion
Wild_df[, "diffusion_low"] <- "km^2/day" %#% Wild_df[, "diffusion_low"] # converts m^2/s to km^2/day
Wild_df[, "diffusion_est"] <- "km^2/day" %#% Wild_df[, "diffusion_est"] # converts m^2/s to km^2/day
Wild_df[, "diffusion_high"] <- "km^2/day" %#% Wild_df[, "diffusion_high"] # converts m^2/s to km^2/day
# home range size estimation 
Wild_df[,"HR_low"] <- "km^2" %#% Wild_df[,"HR_low"] # converts m^2 to km^2
Wild_df[, "HR_est"] <- "km^2" %#% Wild_df[,"HR_est"] # converts m^2 to km^2
Wild_df[, "HR_high"] <- "km^2" %#% Wild_df[,"HR_high"] # converts m^2 to km^2
# mean speed was already in km/day for each individual
#Wild_df[, "Mean_Speed_low"] <- "km/day" %#% Wild_df[, "Mean_Speed_low"] # converts m/s to km/day 
#Wild_df[, "Mean_Speed_est"] <- "km/day" %#% Wild_df[, "Mean_Speed_est"] # converts m/s to km/day 
#Wild_df[, "Mean_Speed_high"] <- "km/day" %#% Wild_df[, "Mean_Speed_high"] # converts m/s to km/day 


# remove rownames
row.names(Wild_df) <- NULL

# add information for each dataframe
Orphan_df$Status <- "Orphaned"
Wild_df$Status <- "Wild-raised"

# add residency infomation (i.e. disperser vs. range-resident)
# orphaned
Orphan_df[-c(2,3,5,6,11,14,15,19,21,22,23,24,26), "Residency"] <- "Range-resident"
Orphan_df[c(2,3,5,6,11,14,15,19,21,22,23,24,26), "Residency"] <- "Disperser"

# wild-raised
Wild_df[-c(24,25,26), "Residency"] <- "Range-resident"
Wild_df[c(24,25,26), "Residency"] <- "Disperser"

# add information about who is successful vs unsuccessful
# orphaned
Orphan_df[c(4,9,10,12,13,15,16), "Success"] <- "Successful"
Orphan_df[c(1,2,5,6,11,14,17,18,19,20,21,22,23,24,25,26,27,28), "Success"] <- "Unsuccessful" 
Orphan_df[c(3,7,8), "Success"] <- "Uncertain" 

#wild-raised
Wild_df[,"Success"] <- NA


# combine the two dataframes into one
Movement_df <- rbind(Orphan_df, Wild_df)

# save file
save(Movement_df , file ="./RESULTS/Movement_model_df.rda")






# meta-analysis ----
# meta function in ctmm calculates means of different movement parameters in a more informed manner 
# (i.e. taking sample size into account) and essentially conducts a meta analysis on the population

# orphaned population
# speed
SPEED <- ctmm::meta(SPEED_orphan, variable = "speed", level = 0.95, units = FALSE)
SPEED <- as.data.frame(SPEED)
SPEED <- SPEED[-c(2,3),] # removes unnecessary rows
SPEED[,"low"] <- "km/day" %#% SPEED[,"low"] # converts m/s to km/day 
SPEED[,"est"] <- "km/day" %#% SPEED[, "est"] # converts m/s to km/day 
SPEED[,"high"] <- "km/day" %#% SPEED[,"high"] # converts m/s to km/day 
colnames(SPEED) <- c("speed_low", "speed_est", "speed_high")

# diffusion
DIFFUSION <- ctmm::meta(FITS_orphan, variable = "diffusion", level = 0.95, units = FALSE)
DIFFUSION <- as.data.frame(DIFFUSION)
DIFFUSION <- DIFFUSION[-c(2,3),] # removes unnecessary rows
DIFFUSION[,"low"] <- "km^2/day" %#% DIFFUSION[,"low"] # converts m^2/s to km^2/day 
DIFFUSION[,"est"] <- "km^2/day" %#% DIFFUSION[, "est"] # converts m^2/s to km^2/day 
DIFFUSION[,"high"] <- "km^2/day" %#% DIFFUSION[,"high"] # converts m^2/s to km^2/day 
colnames(DIFFUSION) <- c("diffusion_low", "diffusion_est", "diffusion_high")

# tau velocity
TAUVELOCITY <- ctmm::meta(FITS_orphan, variable = "tau velocity", level = 0.95, units= FALSE)
TAUVELOCITY <- as.data.frame(TAUVELOCITY)
TAUVELOCITY <- TAUVELOCITY[-c(2,3),] # removes unnecessary rows
TAUVELOCITY[,"low"] <- "minutes" %#% TAUVELOCITY[,"low"] # converts seconds to minutes
TAUVELOCITY[,"est"] <- "minutes" %#% TAUVELOCITY[, "est"] # converts seconds to minutes
TAUVELOCITY[,"high"] <- "minutes" %#% TAUVELOCITY[,"high"] # converts seconds to minutes
colnames(TAUVELOCITY) <- c("τvelocity_low", "τvelocity_est", "τvelocity_high")

# home range size 
HR <- ctmm::meta(AKDE_orphan, variable = "area", level = 0.95, units = FALSE)
HR <- as.data.frame(HR)
HR <- HR[-c(2,3),] # removes unnecessary rows
HR[,"low"] <- "km^2" %#% HR[,"low"] # converts m^2 to km^2
HR[,"est"] <- "km^2" %#% HR[, "est"] # converts m^2 to km^2
HR[,"high"] <- "km^2" %#% HR[,"high"] # converts m^2 to km^2
colnames(HR) <- c("HR_low", "HR_est", "HR_high")
# combine all
META_total_o <- cbind(SPEED, DIFFUSION, TAUVELOCITY, HR)
META_total_o$Status <- "Orphaned"

# wild-raised population
# speed
SPEED <- ctmm::meta(SPEED_wild, variable = "speed", level = 0.95, units = FALSE)
SPEED <- as.data.frame(SPEED)
SPEED <- SPEED[-c(2,3),] # removes unnecessary rows
SPEED[,"low"] <- "km/day" %#% SPEED[,"low"] # converts m/s to km/day 
SPEED[,"est"] <- "km/day" %#% SPEED[, "est"] # converts m/s to km/day 
SPEED[,"high"] <- "km/day" %#% SPEED[,"high"] # converts m/s to km/day 
colnames(SPEED) <- c("speed_low", "speed_est", "speed_high")
# diffusion
DIFFUSION <- ctmm::meta(FITS_wild, variable = "diffusion", level = 0.95, units = FALSE)
DIFFUSION <- as.data.frame(DIFFUSION)
DIFFUSION <- DIFFUSION[-c(2,3),] # removes unnecessary rows
DIFFUSION[,"low"] <- "km^2/day" %#% DIFFUSION[,"low"] # converts m^2/s to km^2/day 
DIFFUSION[,"est"] <- "km^2/day" %#% DIFFUSION[, "est"] # converts m^2/s to km^2/day 
DIFFUSION[,"high"] <- "km^2/day" %#% DIFFUSION[,"high"] # converts m^2/s to km^2/day 
colnames(DIFFUSION) <- c("diffusion_low", "diffusion_est", "diffusion_high")
# tau velocity
TAUVELOCITY <- ctmm::meta(FITS_wild, variable = "tau velocity", level = 0.95, units= FALSE)
TAUVELOCITY <- as.data.frame(TAUVELOCITY)
TAUVELOCITY <- TAUVELOCITY[-c(2,3),] # removes unnecessary rows
TAUVELOCITY[,"low"] <- "minutes" %#% TAUVELOCITY[,"low"] # converts seconds to minutes
TAUVELOCITY[,"est"] <- "minutes" %#% TAUVELOCITY[, "est"] # converts seconds to minutes
TAUVELOCITY[,"high"] <- "minutes" %#% TAUVELOCITY[,"high"] # converts seconds to minutes
colnames(TAUVELOCITY) <- c("τvelocity_low", "τvelocity_est", "τvelocity_high")
# home range
HR <- ctmm::meta(AKDE_wild, variable = "area", level = 0.95, units = FALSE)
HR <- as.data.frame(HR)
HR <- HR[-c(2,3),] # removes unnecessary rows
HR[,"low"] <- "km^2" %#% HR[,"low"] # converts m^2 to km^2
HR[,"est"] <- "km^2" %#% HR[, "est"] # converts m^2 to km^2
HR[,"high"] <- "km^2" %#% HR[,"high"] # converts m^2 to km^2
colnames(HR) <- c("HR_low", "HR_est", "HR_high")
# combine all
META_total_w <- cbind(SPEED, DIFFUSION, TAUVELOCITY, HR)
META_total_w$Status <- "Wild-raised"

# bind orphaned and wild-raised together
META_df <- rbind(META_total_o, META_total_w)

# save output
save(META_df, file = "./RESULTS/META_average_df.rda")






































