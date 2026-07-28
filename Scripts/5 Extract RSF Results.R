# 5 Extract RSF Results
# The main steps carried out by this script are:
#       1. Reorganize RSFs
#       2. Extract result output

#load packages
library(ctmm) #working with ctmm objects
library(dplyr) #working with dataframes

#load RSFs ( several individuals saved yet the loop needed to be restarted)
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2019.rda") #2019 orphaned
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2021.rda") #2019 orphaned
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2022.rda") #2019 orphaned
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2023.rda") #2019 orphaned
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2024.rda") #2019 orphaned
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2017.rda") # wild-raised
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2018.rda") # wild-raised
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2017_Anthony.rda") # wild-raised
Anthony <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2017_Bumpus.rda") # wild-raised
Bumpus <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2017_Cate.rda") # wild-raised
Cate <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2017_Christoffer.rda") # wild-raised
Christoffer <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2017_Elaine.rda") # wild-raised
Elaine <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2017_Jackson.rda") # wild-raised
Jackson <- IND
TEST <- list(Anthony, Bumpus, Cate, Christoffer, Elaine, Jackson)
names(TEST) <- c("Anthony", "Bumpus", "Cate", "Christoffer", "Elaine", "Jackson")
RSF_17 <- c(RSF_17, TEST)

load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2018_Alexander.rda") # wild-raised
Alexander <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2018_Annie.rda") # wild-raised
Annie <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2018_Beto.rda") # wild-raised
Beto <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2018_Hannah.rda") # wild-raised
Hannah <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2018_Jane.rda") # wild-raised
Jane <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2018_Larry.rda") # wild-raised
Larry <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2018_Luigi.rda") # wild-raised
Luigi <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2018_Margaret.rda") # wild-raised
Margaret <- IND
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Individual/2018_Maria.rda") # wild-raised
Maria <- IND
# combine individuals into one list
TEST <- list(Alexander, Annie, Beto, Hannah, Jane, Larry, Luigi, Margaret, Maria)
# add names to list
names(TEST) <- c("Alexander", "Annie", "Beto", "Hannah", "Jane", "Larry", "Luigi", "Margaret", "Maria")
RSF_18 <- c(RSF_18, TEST)
# save results for each year
save(RSF_17, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2017.rda")
save(RSF_18, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_2018.rda")
# combine lists into one for each population
RSF_orphan <- c(RSF_19, RSF_21, RSF_22, RSF_23, RSF_24)
RSF_wild <- c(RSF_17, RSF_18)
#save output
save(RSF_orphan, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_orphan.rda")
save(RSF_wild, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_wild.rda")
#clear up environment space
rm(RSF_19, RSF_21, RSF_22, RSF_23, RSF_17, RSF_18, RSF_24)



# extract individual results -----
# wild-raised -----------------------------------------
#create a dataframe to house total information
RSF_wild_df <- data.frame()
#for.loop for extracting RSF results
for(j in 1:length(RSF_wild)){
  #create a dataframe to house information for each individual
  RSF_extract <- data.frame()
  #extract summary output for an individual
  df <- data.frame(t(summary(RSF_wild[[j]], units = FALSE)[["CI"]]))
  #if value is available in summary output, extract it into the individual's dataframe. 
  #if value is not available in summary output, include an NA in the column
  for(i in 1:1) {
    #Agriculture
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_extract [i, "Agriculture_low"] <- df["low", "Agriculture..1.Agriculture."]
    }else{
      RSF_extract[i, "Agriculture_low"] <- NA
    }
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_extract [i, "Agriculture_est"] <- df["est", "Agriculture..1.Agriculture."]
    }else{
      RSF_extract[i, "Agriculture_est"] <- NA
    }
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_extract [i, "Agriculture_high"] <- df["high", "Agriculture..1.Agriculture."]
    }else{
      RSF_extract[i, "Agriculture_high"] <- NA
    }
    #Development
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_extract [i, "Development_low"] <- df["low", "Development..1.Development."]
    }else{
      RSF_extract[i, "Development_low"] <- NA
    }
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_extract [i, "Development_est"] <- df["est", "Development..1.Development."]
    }else{
      RSF_extract[i, "Development_est"] <- NA
    }
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_extract [i, "Development_high"] <- df["high", "Development..1.Development."]
    }else{
      RSF_extract[i, "Development_high"] <- NA
    }
    #Forest
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_extract [i, "Forest_low"] <- df["low", "Forest..1.Forest."]
    }else{
      RSF_extract[i, "Forest_low"] <- NA
    }
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_extract [i, "Forest_est"] <- df["est", "Forest..1.Forest."]
    }else{
      RSF_extract[i, "Forest_est"] <- NA
    }
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_extract [i, "Forest_high"] <- df["high", "Forest..1.Forest."]
    }else{
      RSF_extract[i, "Forest_high"] <- NA
    }
  
    #Pasture
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_extract [i, "Pasture_low"] <- df["low", "Pasture..1.Pasture."]
    }else{
      RSF_extract[i, "Pasture_low"] <- NA
    }
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_extract [i, "Pasture_est"] <- df["est", "Pasture..1.Pasture."]
    }else{
      RSF_extract[i, "Pasture_est"] <- NA
    }
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_extract [i, "Pasture_high"] <- df["high", "Pasture..1.Pasture."]
    }else{
      RSF_extract[i, "Pasture_high"] <- NA
    }
    #Water
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_extract [i, "Water_low"] <- df["low", "Water..1.Water."]
    }else{
      RSF_extract[i, "Water_low"] <- NA
    }
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_extract [i, "Water_est"] <- df["est", "Water..1.Water."]
    }else{
      RSF_extract[i, "Water_est"] <- NA
    }
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_extract [i, "Water_high"] <- df["high", "Water..1.Water."]
    }else{
      RSF_extract[i, "Water_high"] <- NA
    }} #end of inner loop
  #add ID of individual
  RSF_extract$ID <- RSF_wild[[j]]@info$identity
  #bind individual results to total dataframe
  RSF_wild_df <- rbind(RSF_wild_df, RSF_extract)
  
} #end of outer loop

#save output
save(RSF_wild_df, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_wild_df.rda")




# orphaned results
#create a dataframe to hold total results
RSF_orphan_df <- data.frame()
#for.loop for extracting RSF results
for(j in 1:length(RSF_orphan)){
  #create a dataframe to hold an individual's results
  RSF_extract <- data.frame()
  #extract summary output for an individual
  df <- data.frame(t(summary(RSF_orphan[[j]], units = FALSE)[["CI"]]))
  #if value is available in summary output, extract it into the individual's dataframe. 
  #if value is not available in summary output, include an NA in the column
  for(i in 1:1) {
    #Agriculture
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_extract [i, "Agriculture_low"] <- df["low", "Agriculture..1.Agriculture."]
    }else{
      RSF_extract[i, "Agriculture_low"] <- NA
    }
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_extract [i, "Agriculture_est"] <- df["est", "Agriculture..1.Agriculture."]
    }else{
      RSF_extract[i, "Agriculture_est"] <- NA
    }
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_extract [i, "Agriculture_high"] <- df["high", "Agriculture..1.Agriculture."]
    }else{
      RSF_extract[i, "Agriculture_high"] <- NA
    }
    #Development
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_extract [i, "Development_low"] <- df["low", "Development..1.Development."]
    }else{
      RSF_extract[i, "Development_low"] <- NA
    }
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_extract [i, "Development_est"] <- df["est", "Development..1.Development."]
    }else{
      RSF_extract[i, "Development_est"] <- NA
    }
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_extract [i, "Development_high"] <- df["high", "Development..1.Development."]
    }else{
      RSF_extract[i, "Development_high"] <- NA
    }
    #Forest
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_extract [i, "Forest_low"] <- df["low", "Forest..1.Forest."]
    }else{
      RSF_extract[i, "Forest_low"] <- NA
    }
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_extract [i, "Forest_est"] <- df["est", "Forest..1.Forest."]
    }else{
      RSF_extract[i, "Forest_est"] <- NA
    }
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_extract [i, "Forest_high"] <- df["high", "Forest..1.Forest."]
    }else{
      RSF_extract[i, "Forest_high"] <- NA
    }
    
    #Pasture
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_extract [i, "Pasture_low"] <- df["low", "Pasture..1.Pasture."]
    }else{
      RSF_extract[i, "Pasture_low"] <- NA
    }
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_extract [i, "Pasture_est"] <- df["est", "Pasture..1.Pasture."]
    }else{
      RSF_extract[i, "Pasture_est"] <- NA
    }
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_extract [i, "Pasture_high"] <- df["high", "Pasture..1.Pasture."]
    }else{
      RSF_extract[i, "Pasture_high"] <- NA
    }
    #Water
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_extract [i, "Water_low"] <- df["low", "Water..1.Water."]
    }else{
      RSF_extract[i, "Water_low"] <- NA
    }
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_extract [i, "Water_est"] <- df["est", "Water..1.Water."]
    }else{
      RSF_extract[i, "Water_est"] <- NA
    }
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_extract [i, "Water_high"] <- df["high", "Water..1.Water."]
    }else{
      RSF_extract[i, "Water_high"] <- NA
    }} #end of inner loop
  
  #add individual's ID
  RSF_extract$ID <- RSF_orphan[[j]]@info$identity
  #add individual dataframe to total dataframe
  RSF_orphan_df <- rbind(RSF_orphan_df, RSF_extract) 
} #end of outer loop

#save output
save(RSF_orphan_df, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_orphan_df.rda")




# calculate population averages for selection ------
RSF_meta_wild <- mean(RSF_wild)
RSF_meta_orphan <- mean(RSF_orphan)
# save results
save(RSF_mean_orphan, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Mean_RSF_orphan.rda")
save(RSF_mean_wild, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Mean_RSF_wild.rda")

# wild-raised population
# create a dataframe to hold an individual's results
  RSF_mean_wild_df <- data.frame()
  
  # extract summary output for an individual
  df <- data.frame(t(summary(RSF_mean_wild, units = FALSE)[["CI"]]))
  # if value is available in summary output, extract it into the individual's dataframe. 
  # if value is not available in summary output, include an NA in the column
  for(i in 1:1) {
    #Agriculture
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Agriculture_low"] <- df["low", "Agriculture..1.Agriculture."]
    }else{
      RSF_mean_wild_df[i, "Agriculture_low"] <- NA
    }
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Agriculture_est"] <- df["est", "Agriculture..1.Agriculture."]
    }else{
      RSF_mean_wild_df[i, "Agriculture_est"] <- NA
    }
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Agriculture_high"] <- df["high", "Agriculture..1.Agriculture."]
    }else{
      RSF_mean_wild_df[i, "Agriculture_high"] <- NA
    }
    #Development
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Development_low"] <- df["low", "Development..1.Development."]
    }else{
      RSF_mean_wild_df[i, "Development_low"] <- NA
    }
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Development_est"] <- df["est", "Development..1.Development."]
    }else{
      RSF_mean_wild_df[i, "Development_est"] <- NA
    }
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Development_high"] <- df["high", "Development..1.Development."]
    }else{
      RSF_mean_wild_df[i, "Development_high"] <- NA
    }
    #Forest
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Forest_low"] <- df["low", "Forest..1.Forest."]
    }else{
      RSF_mean_wild_df[i, "Forest_low"] <- NA
    }
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Forest_est"] <- df["est", "Forest..1.Forest."]
    }else{
      RSF_mean_wild_df[i, "Forest_est"] <- NA
    }
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Forest_high"] <- df["high", "Forest..1.Forest."]
    }else{
      RSF_mean_wild_df[i, "Forst_high"] <- NA
    }
    
    #Pasture
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Pasture_low"] <- df["low", "Pasture..1.Pasture."]
    }else{
      RSF_mean_wild_df[i, "Pasture_low"] <- NA
    }
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Pasture_est"] <- df["est", "Pasture..1.Pasture."]
    }else{
      RSF_mean_wild_df[i, "Pasture_est"] <- NA
    }
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Pasture_high"] <- df["high", "Pasture..1.Pasture."]
    }else{
      RSF_mean_wild_df[i, "Pasture_high"] <- NA
    }
    #Water
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Water_low"] <- df["low", "Water..1.Water."]
    }else{
      RSF_mean_wild_df[i, "Water_low"] <- NA
    }
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Water_est"] <- df["est", "Water..1.Water."]
    }else{
      RSF_mean_wild_df[i, "Water_est"] <- NA
    }
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_mean_wild_df [i, "Water_high"] <- df["high", "Water..1.Water."]
    }else{
      RSF_mean_wild_df[i, "Water_high"] <- NA
    }} #end of inner loop
  
  #add individual's ID
  RSF_mean_wild_df$ID <- "Wild-raised"
  
  save(RSF_mean_wild_df, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Mean_RSF_wild_df.rda")
  

  
  # orphaned population ----
  # create a dataframe to hold an individual's results
  RSF_mean_orphan_df <- data.frame()
  
  # extract summary output for an individual
  df <- data.frame(t(summary(RSF_mean_orphan, units = FALSE)[["CI"]]))
  # if value is available in summary output, extract it into the individual's dataframe. 
  # if value is not available in summary output, include an NA in the column
  for(i in 1:1) {
    #Agriculture
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Agriculture_low"] <- df["low", "Agriculture..1.Agriculture."]
    }else{
      RSF_mean_orphan_df[i, "Agriculture_low"] <- NA
    }
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Agriculture_est"] <- df["est", "Agriculture..1.Agriculture."]
    }else{
      RSF_mean_orphan_df[i, "Agriculture_est"] <- NA
    }
    if ("Agriculture..1.Agriculture." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Agriculture_high"] <- df["high", "Agriculture..1.Agriculture."]
    }else{
      RSF_mean_orphan_df[i, "Agriculture_high"] <- NA
    }
    #Development
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Development_low"] <- df["low", "Development..1.Development."]
    }else{
      RSF_mean_orphan_df[i, "Development_low"] <- NA
    }
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Development_est"] <- df["est", "Development..1.Development."]
    }else{
      RSF_mean_orphan_df[i, "Development_est"] <- NA
    }
    if ("Development..1.Development." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Development_high"] <- df["high", "Development..1.Development."]
    }else{
      RSF_mean_orphan_df[i, "Development_high"] <- NA
    }
    #Forest
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Forest_low"] <- df["low", "Forest..1.Forest."]
    }else{
      RSF_mean_orphan_df[i, "Forest_low"] <- NA
    }
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Forest_est"] <- df["est", "Forest..1.Forest."]
    }else{
      RSF_mean_orphan_df[i, "Forest_est"] <- NA
    }
    if ("Forest..1.Forest." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Forest_high"] <- df["high", "Forest..1.Forest."]
    }else{
      RSF_mean_orphan_df[i, "Forst_high"] <- NA
    }
    
    #Pasture
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Pasture_low"] <- df["low", "Pasture..1.Pasture."]
    }else{
      RSF_mean_orphan_df[i, "Pasture_low"] <- NA
    }
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Pasture_est"] <- df["est", "Pasture..1.Pasture."]
    }else{
      RSF_mean_orphan_df[i, "Pasture_est"] <- NA
    }
    if ("Pasture..1.Pasture." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Pasture_high"] <- df["high", "Pasture..1.Pasture."]
    }else{
      RSF_mean_orphan_df[i, "Pasture_high"] <- NA
    }
    #Water
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Water_low"] <- df["low", "Water..1.Water."]
    }else{
      RSF_mean_orphan_df[i, "Water_low"] <- NA
    }
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Water_est"] <- df["est", "Water..1.Water."]
    }else{
      RSF_mean_orphan_df[i, "Water_est"] <- NA
    }
    if ("Water..1.Water." %in% colnames(df)) {
      RSF_mean_orphan_df [i, "Water_high"] <- df["high", "Water..1.Water."]
    }else{
      RSF_mean_orphan_df[i, "Water_high"] <- NA
    }} #end of inner loop
  
  # add individual's ID
  RSF_mean_orphan_df$ID <- "Orphaned"
  # save results
  save(RSF_mean_orphan_df, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Mean_RSF_wild_df.rda")
  
  # combine results into one dataframe
  mean_all <- rbind(RSF_mean_orphan_df, RSF_mean_wild_df)  
  # save results
  save(mean_all, file = "./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Mean_RSF_all_df.rda")
