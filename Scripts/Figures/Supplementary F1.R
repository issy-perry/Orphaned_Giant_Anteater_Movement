# load packages
library(ctmm) #movement
library(dplyr)
library(ggplot2) #making plots
library(ggpubr) #arranging plots in multiples
library(tidyr) #drop_na function
library(patchwork) # plotting

# load data
# UDs
load("./RESULTS/AKDEs/UDs_wild_RR.rda")
# model fits
load("./RESULTS/Fits/Fits_wild_RR.rda")
# Speeds
load("./RESULTS/Speed/Speed_wild_raised_RR.rda")

# subadults are anthony, gala, kyle, little_rick, and reid
SUBADULT_FIT <- FITS_wild [c("Anthony", "Kyle", "Little_Rick", "Reid")]
SUBADULT_AKDE <- AKDE_wild[c("Anthony", "Kyle", "Little_Rick", "Reid")]

ADULT_FIT <- FITS_wild[c("Bumpus", "Cate", "Christoffer", "Elaine", "Jackson", "Makao", 
                         "Puji", "Alexander", "Annie", "Beto","Hannah", "Jane", "Larry", 
                         "Luigi", "Margaret", "Maria","Rodolfo", "Sheron", "Thomas")]
ADULT_AKDE <- AKDE_wild[c("Bumpus", "Cate", "Christoffer", "Elaine", "Jackson", "Makao", 
                          "Puji", "Alexander", "Annie", "Beto","Hannah", "Jane", "Larry", 
                          "Luigi", "Margaret", "Maria","Rodolfo", "Sheron", "Thomas")]



# load in total dataframe for plotting
load("RESULTS/Movement_model_RR_df.rda")
# add age information to dataframe
Movement_df[Movement_df$ID %in% c("Bumpus", "Cate", "Christoffer", "Elaine", "Jackson", "Makao", 
              "Puji", "Alexander", "Annie", "Beto","Hannah", "Jane", "Larry", 
              "Luigi", "Margaret", "Maria","Rodolfo", "Sheron", "Thomas"), "Status"] <- "Adult"

Movement_df[Movement_df$ID %in% c("Anthony", "Kyle", "Little_Rick", "Reid"), "Status"] <- "Subadult"

# total averages ----

# calculate movement metric averages for each age class
#wild subadult
HR_J <- ctmm::meta(SUBADULT_AKDE, variable = "area", level = 0.95, units = FALSE)
HR_J <- as.data.frame(HR_J)
HR_J <- HR_J[-c(2,3),]
HR_J[,"low"] <- "km^2" %#% HR_J[,"low"] # converts to km/day rather than m/s
HR_J[,"est"] <- "km^2" %#% HR_J[, "est"] # converts to km/day rather than m/s
HR_J[,"high"] <- "km^2" %#% HR_J[,"high"] # converts to km/day rather than m/s
colnames(HR_J) <- c("HR_low", "HR_est", "HR_high")

#wild adult
HR_A <- ctmm::meta(ADULT_AKDE, variable = "area", level = 0.95, units = FALSE)
HR_A <- as.data.frame(HR_A)
HR_A <- HR_A[-c(2,3),]
HR_A[,"low"] <- "km^2" %#% HR_A[,"low"] # converts to km/day rather than m/s
HR_A[,"est"] <- "km^2" %#% HR_A[, "est"] # converts to km/day rather than m/s
HR_A[,"high"] <- "km^2" %#% HR_A[,"high"] # converts to km/day rather than m/s
colnames(HR_A) <- c("HR_low", "HR_est", "HR_high")



#wild subadult
SPEED_J <- ctmm::meta(SUBADULT_FIT, variable = "speed", level = 0.95, units = FALSE)
SPEED_J <- as.data.frame(SPEED_J)
SPEED_J <- SPEED_J[-c(2,3),]
SPEED_J[,"low"] <- "km/day" %#% SPEED_J[,"low"] # converts to km/day rather than m/s
SPEED_J[,"est"] <- "km/day" %#% SPEED_J[, "est"] # converts to km/day rather than m/s
SPEED_J[,"high"] <- "km/day" %#% SPEED_J[,"high"] # converts to km/day rather than m/s
colnames(SPEED_J) <- c("speed_low", "speed_est", "speed_high")

#wild adult
SPEED_A <- ctmm::meta(ADULT_FIT, variable = "speed", level = 0.95, units = FALSE)
SPEED_A <- as.data.frame(SPEED_A)
SPEED_A <- SPEED_A[-c(2,3),]
SPEED_A[,"low"] <- "km/day" %#% SPEED_A[,"low"] # converts to km/day rather than m/s
SPEED_A[,"est"] <- "km/day" %#% SPEED_A[, "est"] # converts to km/day rather than m/s
SPEED_A[,"high"] <- "km/day" %#% SPEED_A[,"high"] # converts to km/day rather than m/s
colnames(SPEED_A) <- c("speed_low", "speed_est", "speed_high")


#wild subadult
VELOCITY_J <- ctmm::meta(SUBADULT_FIT, variable = "tauvelocity", level = 0.95, units = FALSE)
VELOCITY_J <- as.data.frame(VELOCITY_J)
VELOCITY_J <- VELOCITY_J[-c(2,3),]
VELOCITY_J[,"low"] <- "minutes" %#% VELOCITY_J[,"low"] # converts to km/day rather than m/s
VELOCITY_J[,"est"] <- "minutes" %#% VELOCITY_J[, "est"] # converts to km/day rather than m/s
VELOCITY_J[,"high"] <- "minutes" %#% VELOCITY_J[,"high"] # converts to km/day rather than m/s
colnames(VELOCITY_J) <- c("τvelocity_low", "τvelocity_est", "τvelocity_high")

#wild adult
VELOCITY_A <- ctmm::meta(ADULT_FIT, variable = "tauvelocity", level = 0.95, units = FALSE)
VELOCITY_A <- as.data.frame(VELOCITY_A)
VELOCITY_A <- VELOCITY_A[-c(2,3),]
VELOCITY_A[,"low"] <- "minutes" %#% VELOCITY_A[,"low"] # converts to km/day rather than m/s
VELOCITY_A[,"est"] <- "minutes" %#% VELOCITY_A[, "est"] # converts to km/day rather than m/s
VELOCITY_A[,"high"] <- "minutes" %#% VELOCITY_A[,"high"] # converts to km/day rather than m/s
colnames(VELOCITY_A) <- c("τvelocity_low", "τvelocity_est", "τvelocity_high")

# bind results for each age class together
JUV <- cbind(HR_J, SPEED_J, VELOCITY_J)
ADULT <- cbind(HR_A, SPEED_A, VELOCITY_A)
# add informatin about age class
JUV$Status <- "Subadult"
ADULT$Status <- "Adult"
# bind both age class results together
META_df <- rbind(JUV, ADULT)
# free up environment space
rm(HR_J, SPEED_J, VELOCITY_J, HR_A, SPEED_A, VELOCITY_A, WILD)

# reorder IDs in alpbetical order
Movement_df <- Movement_df %>%
  arrange(ID)

# convert ID values to numeric for plotting (needs to be numeric when on the y-axis)
Movement_df$ID <- as.character(Movement_df$ID) 

# plot differences in home range size 
HR <- ggplot() +
  # adult
  geom_pointrange(data = subset(Movement_df, Status == "Adult"), 
                  aes(x = HR_est, y = ID, xmin = HR_low, xmax = HR_high, color = "Adult"), 
                  shape = 17) +
  geom_vline(xintercept = META_df[2,"HR_est"],
             color = "#8732D2", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[2,], 
            aes(xmin = HR_low, xmax = HR_high, ymin = -Inf, ymax = Inf, color = "Adult"), 
            color = NA, fill = "#8732D2", alpha = 0.15) +
  # subadult
  geom_pointrange(data = subset(Movement_df, Status == "Subadult"), 
                  aes(x = HR_est, y = ID, xmin = HR_low, xmax = HR_high, color = "Subadult"), 
                  shape = 15) +
  geom_vline(xintercept = META_df[1,"HR_est"], color = "#DA1C1C", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[1,], 
            aes(xmin = HR_low, xmax = HR_high, ymin = -Inf, ymax = Inf, color = "Subadult"), 
            color = NA, fill = "#DA1C1C", alpha = 0.1) +
  scale_color_manual("Age Class", 
                     values = c("Subadult" = "#DA1C1C", "Adult" = "#8732D2")) +
  labs(title = "a.", 
       y = "", 
       x = "log(Home Range (km^2))") +
  scale_x_continuous(trans = "log10") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "darkgray"), 
        plot.title = element_text(size = 12, family = "sans", hjust = 0.05),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 12),
        legend.position = "none",
        legend.text = element_text(size=12, family = "sans", face = "bold"),)  
#plot(HR)

# plot differences in directional persistence
tauvelocity <- ggplot() +
  # adult
  geom_pointrange(data = subset(Movement_df, Status == "Adult"), 
                  aes(x = τvelocity_est, y = ID, xmin = τvelocity_low, xmax = τvelocity_high, color = "Adult"), 
                  shape = 17) +
  geom_vline(xintercept = META_df[2,"τvelocity_est"],
             color = "#8732D2", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[2,], 
            aes(xmin = τvelocity_low, xmax = τvelocity_high, ymin = -Inf, ymax = Inf, color = "Adult"), 
            color = NA, fill = "#8732D2", alpha = 0.15) +
  # subadult
  geom_pointrange(data = subset(Movement_df, Status == "Subadult"), 
                  aes(x = τvelocity_est, y = ID, xmin = τvelocity_low, xmax = τvelocity_high, color = "Subadult"), 
                  shape = 15) +
  geom_vline(xintercept = META_df[1,"τvelocity_est"], color = "#DA1C1C", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[1,], 
            aes(xmin = τvelocity_low, xmax = τvelocity_high, ymin = -Inf, ymax = Inf, color = "Subadult"), 
            color = NA, fill = "#DA1C1C", alpha = 0.1) +
  scale_color_manual("Age Class", 
                     values = c("Subadult" = "#DA1C1C", "Adult" = "#8732D2")) +
  labs(title = "b.", 
       y = "", 
       x = "τvelocity (minutes)") +
  #scale_x_continuous(trans = "log10") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "darkgray"), 
        plot.title = element_text(size = 12, family = "sans", hjust = 0.05),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 12),
        legend.position = "none",
        legend.text = element_text(size=12, family = "sans", face = "bold"),)  
#plot(tauvelocity)


# plot differences in speed
speed <- ggplot() +
  # adult
  geom_pointrange(data = subset(Movement_df, Status == "Adult" & is.finite(Mean_Speed_est) & is.finite(Mean_Speed_low) & is.finite(Mean_Speed_high)), 
                  aes(x = Mean_Speed_est, y = ID, xmin = Mean_Speed_low, xmax = Mean_Speed_high, color = "Adult"), 
                  shape = 17) +
  geom_vline(xintercept = META_df[2,"speed_est"],
             color = "#8732D2", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[2,], 
            aes(xmin = speed_low, xmax = speed_high, ymin = -Inf, ymax = Inf, color = "Adult"), 
            color = NA, fill = "#8732D2", alpha = 0.15) +
  # subadult
  geom_pointrange(data = subset(Movement_df, Status == "Subadult" & is.finite(Mean_Speed_est) & is.finite(Mean_Speed_low) & is.finite(Mean_Speed_high)), 
                  aes(x = Mean_Speed_est, y = ID, xmin = Mean_Speed_low, xmax = Mean_Speed_high, color = "Subadult"), 
                  shape = 15) +
  geom_vline(xintercept = META_df[1,"speed_est"], color = "#DA1C1C", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[1,], 
            aes(xmin = speed_low, xmax = speed_high, ymin = -Inf, ymax = Inf, color = "Subadult"), 
            color = NA, fill = "#DA1C1C", alpha = 0.1) +
  scale_color_manual("Age Class", 
                     values = c("Subadult" = "#DA1C1C", "Adult" = "#8732D2")) +
  labs(title = "c.", 
       y = "", 
       x = "Speed (km/day)") +
  #scale_x_continuous(trans = "log10") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "darkgray"), 
        plot.title = element_text(size = 12, family = "sans", hjust = 0.05),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=12, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 12),)  
#plot(speed)



# window ratio ----

# load telemetry data for names
load("./DATA/Wild_raised_tel_data/Data_telemetry.rda") 

# separate based on age
DATA_subadults <- DATA_wild[c("Anthony", "Kyle", "Little_Rick", "Reid")]
DATA_adults <- DATA_wild[c("Alexander", "Annie", "Beto", "Bumpus", "Cate", "Christoffer", "Elaine", "Hannah", "Jackson", "Jane",
                      "Larry", "Luigi", "Makao", "Margaret", "Maria", "Puji", "Rodolfo", "Sheron", "Thomas")]


# load comparison data of population average
load("./RESULTS/Fits/Fits_wild.rda")
load("./RESULTS/Speed/Speed_wild_raised.rda")
# separate based on age class
FITS_adults <- FITS_wild[c("Alexander", "Annie", "Beto", "Bumpus", "Cate", "Christoffer", "Elaine", "Hannah", "Jackson", "Jane",
                           "Larry", "Luigi", "Makao", "Margaret", "Maria", "Puji", "Rodolfo", "Sheron", "Thomas")]
SPEED_adults <- SPEED_wild[c("Alexander", "Annie", "Beto", "Bumpus", "Cate", "Christoffer", "Elaine", "Hannah", "Jackson", "Jane",
                            "Larry", "Luigi", "Makao", "Margaret", "Maria", "Puji", "Rodolfo", "Sheron", "Thomas")]


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
fits_s <- file.path(fits_folder, paste0("FIT_", names(DATA_subadults), ".rda")) 

# pull fits for orphans
FITS_subadults <- lapply(fits_s, function(cur_file) {
  load(cur_file)
  return(FITS)
})

# rename the lists based on identity rather than another named list as these are now alphabetical
for(i in 1:length(FITS_subadults)) {
  FITS_ID <- FITS_subadults[[i]]
  names(FITS_subadults)[i] <- FITS_ID[[1]]@info$identity
}


# AKDEs
akde_s <- file.path(akde_folder, paste0("UD_", names(DATA_subadults), ".rda")) 

# pull UDs for orphans
AKDE_subadults <- lapply(akde_s, function(cur_file) {
  load(cur_file)
  return(wAKDEs)
})

# rename the lists based on identity rather than another named list as these are now alphabetical
for(i in 1:length(AKDE_subadults)) {
  FITS_ID <- FITS_subadults[[i]]
  names(AKDE_subadults)[i] <- FITS_ID[[1]]@info$identity
}


# make file path of wild-raised individuals
akde_w <- file.path(akde_folder, paste0("UD_", names(DATA_wild), ".rda"))

# load UDs for wild-raised individuals
AKDE_wild <- lapply(akde_w, function(cur_file) {
  load(cur_file)
  return(wAKDEs)
})


# Speeds
# make file path of orphans
speed_s <- file.path(speed_folder, paste0("Speed_", names(DATA_subadults), ".rda")) 

# pull speeds for orphans
SPEED_subadults <- lapply(speed_s, function(cur_file) {
  load(cur_file)
  return(speed_mean)
})

# rename the lists based on identity rather than another named list as these are now alphabetical
for(i in 1:length(SPEED_subadults)) {
  FITS_ID <- FITS_subadults[[i]]
  names(SPEED_subadults)[i] <- FITS_ID[[1]]@info$identity
}


# for.loop for extracting age class averages as monitoring progresses
# create a dataframe to house the information
META_win_subadults <- data.frame()

# for.loop pulls windows from individuals
for(j in 1:167){ # number chosen as it is the third highest number of lists in the nested list
  
  # make list to hold extracted windows
  fits_win <- list()
  
  # work within a new for.loop to pull windows from each individual 
  for(i in 1:length(FITS_subadults)){
    
    # pull out individual 
    IND <- FITS_subadults[[i]] 
    
    # pull out window and stick it into a slot on new list
    if(j <= length(IND)) {fits_win[[length(fits_win)+1]] <- IND[[j]]} 
  } # end inner inner loop
  
  #extract AKDEs
  akde_win <- list()
  for(i in 1:length(AKDE_subadults)){
    # pull out individual 
    IND <- AKDE_subadults[[i]] 
    
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
  for(i in 1:length(SPEED_subadults)){
    # pull out individual 
    IND <- SPEED_subadults[[i]]
    
    # pull out window and stick it into a slot on new list
    if(j <= length(IND)) {speed_win[[length(speed_win)+1]] <- IND[[j]]}
  } # end inner inner loop
  
  # ensure only non fractal models are included in analysis
  akde_win <- akde_win[sapply(akde_win, function(x) length(x) >= 13)] # non fractal wAKDE object has 13 elements
  akde_wild <- akde_wild[sapply(akde_wild, function(x) length(x) >= 13)] # non fractal wAKDE object has 13 elements
  speed_win <- speed_win[sapply(speed_win, function(x) {x[[1]] != 0})]  # DOF values should not be 0
  
  
  # home range
  if (length(akde_win) > 2) {
    HR <- ctmm::meta(list(subadult = akde_win, adult = akde_wild), variable = "area", level = 0.95, units = FALSE) # units don't matter, it's a ratio
    HR_o <- as.data.frame(HR)
    HR_o$HR_low <- HR_o["subadult/", "/adult.low"]
    HR_o$HR_est <- HR_o["subadult/", "/adult.est"]
    HR_o$HR_high <- HR_o["subadult/", "/adult.high"]
    HR_o <- HR_o[-c(2), -c(1,2,3,4,5,6)] # removes column 2 and first six rows
  } else {
    HR_o <- data.frame(HR_low = NA,
                       HR_est = NA,
                       HR_high = NA)
  }
  
  # diffusion
  if (length(fits_win) > 2) {
    DIFFUSION <- ctmm::meta(list(subadult = fits_win, adult = FITS_wild), variable = "diffusion", level = 0.95, units = FALSE) # units don't matter, it's a ratio
    DIFFUSION_o <- as.data.frame(DIFFUSION)
    DIFFUSION_o$diffusion_low <- DIFFUSION_o["subadult/", "/adult.low"]
    DIFFUSION_o$diffusion_est <- DIFFUSION_o["subadult/", "/adult.est"]
    DIFFUSION_o$diffusion_high <- DIFFUSION_o["subadult/", "/adult.high"]
    DIFFUSION_o <- DIFFUSION_o[-c(2), -c(1,2,3,4,5,6)] # removes row two and first six columns 
  } else {
    DIFFUSION_o <- data.frame(diffusion_low = NA,
                              diffusion_est = NA,
                              diffusion_high= NA)
  }
  # speed
  if (length(speed_win) > 2) {
    SPEED <- ctmm::meta(list(subadult = speed_win, adult = SPEED_wild), variable = "speed", level = 0.95, units = FALSE) # units don't matter, it's a ratio
    SPEED_o <- as.data.frame(SPEED)
    SPEED_o$speed_low <- SPEED_o["subadult/", "/adult.low"]
    SPEED_o$speed_est <- SPEED_o["subadult/", "/adult.est"]
    SPEED_o$speed_high <- SPEED_o["subadult/", "/adult.high"]
    SPEED_o <- SPEED_o[-c(2), -c(1,2,3,4,5,6)] # removes row two and first six columns 
  } else {
    SPEED_o <- data.frame(speed_low = NA,
                          speed_est = NA,
                          speed_high = NA)
  }
  # tau velocity
  tryCatch({
    if (length(fits_win) > 2) {
      TAUV <- ctmm::meta(list(subadult = fits_win, adult = FITS_wild), variable = "tau velocity", level = 0.95, units = FALSE) # units don't matter, it's a ratio
      TAUV_o <- as.data.frame(TAUV)
      TAUV_o$tauvelocity_low <- TAUV_o["subadult/", "/adult.low"]
      TAUV_o$tauvelocity_est <- TAUV_o["subadult/", "/adult.est"]
      TAUV_o$tauvelocity_high <- TAUV_o["subadult/", "/adult.high"]
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
  META_win_subadults <- rbind(META_win_subadults, META_o) 
  
}#end of outer loop

# remove rownames
rownames(META_win_subadults) <- NULL

# remove rows with NA values for all values
META_win_subadults <- META_win_subadults[!is.na(META_win_subadults$HR_low),]

# add information regarding dNULL# add information regarding days since release
META_win_subadults$time_since_release = seq(10, by = 2, length.out = nrow(META_win_subadults)) # days since release for plotting 

# save output
save(META_win_subadults, file = "./RESULTS/Window/Supplementary_meta_subadults.rda")



HR_win <- ggplot() +
  geom_line(data = META_win_subadults, 
            aes(x = time_since_release, y = HR_est, color = "Ratio of Subadults to Adults"), 
            linewidth = 1) +
  geom_ribbon(data = META_win_subadults, 
              aes (x = time_since_release, ymin = HR_low, ymax = HR_high, fill = "Ratio of Subadults to Adults"), 
              alpha = 0.1) +
  geom_hline(yintercept = 1, 
             col = "grey", 
             linetype = "dashed",
             size = 0.75) +
  #scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c("Ratio of Subadults to Adults" = "#DA1C1C")) + 
  scale_color_manual(values = c("Ratio of Subadults to Adults" = "#DA1C1C")) + 
  labs(title = "d.", 
       y = "home range[subadult] / home range[adult]", 
       x = "Time Since Release (days)") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "darkgray"), 
        legend.position = "none",
        plot.title = element_text(size = 12, family = "sans", hjust = 0.06),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10)) #remove if you want the grid
#plot(HR_win)


speed_win <- ggplot() +
  geom_line(data = META_win_subadults, 
            aes(x = time_since_release, y = speed_est, color = "Ratio of Subadults to Adults"), 
            linewidth = 1) +
  geom_ribbon(data = META_win_subadults, 
              aes (x = time_since_release, ymin = speed_low, ymax = speed_high, fill = "Ratio of Subadults to Adults"), 
              alpha = 0.1) +
  geom_hline(yintercept = 1, 
             col = "grey", 
             linetype = "dashed",
             size = 0.75) +
  #scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c("Ratio of Subadults to Adults" = "#DA1C1C")) + 
  scale_color_manual(values = c("Ratio of Subadults to Adults" = "#DA1C1C")) + 
  labs(title = "f.", y = "speed[subadult] / speed[adult]",
       x = "Time Since Release (days)") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "darkgray"), 
        legend.position = "none",
        plot.title = element_text(size = 12, family = "sans", hjust = 0.06),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10)) #remove if you want the grid
#plot(speed_win)



tauv_win <- ggplot() +
  geom_line(data = META_win_subadults, 
            aes(x = time_since_release, y = tauvelocity_est, color = "Ratio of Subadults to Adults"), 
            linewidth = 1) +
  geom_ribbon(data = META_win_subadults, 
              aes (x = time_since_release, ymin = tauvelocity_low, ymax = tauvelocity_high, fill = "Ratio of Subadults to Adults"), 
              alpha = 0.1) +
  geom_hline(yintercept = 1, 
             col = "grey", 
             linetype = "dashed",
             size = 0.75) +
  #scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c("Ratio of Subadults to Adults" = "#DA1C1C")) + 
  scale_color_manual(values = c("Ratio of Subadults to Adults" = "#DA1C1C")) + 
  labs(title = "e.", 
       y = "tau velocity[subadult] / tauvelocity[adult]", 
       x = "Time Since Release (days)") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "darkgray"), 
        legend.position = "none",
        plot.title = element_text(size = 12, family = "sans", hjust = 0.06),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10)) #remove if you want the grid
#plot(tauv_win)

# create layout for combination of plots
layout_design <- "
  ABC
  DEF
"

# combine all plots together
HR + tauvelocity + speed +
  HR_win + tauv_win + speed_win + plot_layout(design = layout_design) &   
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.title.position = "plot"
  )

# save final plot
ggsave("./FIGURES/Revised/Figure_S1.png", width = 8.5, height = 8, units = "in", 
       dpi = 300, bg = "transparent")

