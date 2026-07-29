# load packages
library(ctmm)
library(ggplot2)
library(ggpubr)
library(dplyr)


# results
load("./RESULTS/META_average_RR_df.rda")
load("./RESULTS/Movement_model_RR_df.rda")
load("./RESULTS/Fits/Fits_orphan_RR.rda")
load("./RESULTS/AKDEs/UDs_orphan_RR.rda") 
load("./RESULTS/Speed/Speed_orphans_RR.rda") 


# separate results based on sex
FO_AKDE <- AKDE_orphan[c("Arya", "Dumbo_1", "Dumbo_2", "Colete", 
                         "Heather", "Juju_2", "Mulan", "Rita", "Bella", "Nancy")]
FO_FITS <- FITS_orphan[c("Arya", "Dumbo_1", "Dumbo_2", "Colete", 
                         "Heather", "Juju_2", "Mulan", "Rita", "Bella", "Nancy")]
FO_SPEED <- SPEED_orphan[c("Arya", "Dumbo_1", "Dumbo_2", "Colete", 
                           "Heather", "Juju_2", "Mulan", "Rita", "Bella", "Nancy")]
MO_AKDE <- AKDE_orphan[c("Peter", "George", "Dom", "Tim_3")]
MO_FITS <- FITS_orphan[c("Peter", "George", "Dom", "Tim_3")]
MO_SPEED <- SPEED_orphan[c("Peter", "George", "Dom", "Tim_3")]

# get averages for each sex
# female orphaned populations
# speed
SPEED <- ctmm::meta(FO_SPEED, variable = "speed", level = 0.95, units = FALSE)
SPEED <- as.data.frame(SPEED)
SPEED <- SPEED[-c(2,3),] # removes unnecessary rows
SPEED[,"low"] <- "km/day" %#% SPEED[,"low"] # converts m/s to km/day 
SPEED[,"est"] <- "km/day" %#% SPEED[, "est"] # converts m/s to km/day 
SPEED[,"high"] <- "km/day" %#% SPEED[,"high"] # converts m/s to km/day 
colnames(SPEED) <- c("speed_low", "speed_est", "speed_high")

# home range size 
HR <- ctmm::meta(FO_AKDE, variable = "area", level = 0.95, units = FALSE)
HR <- as.data.frame(HR)
HR <- HR[-c(2,3),] # removes unnecessary rows
HR[,"low"] <- "km^2" %#% HR[,"low"] # converts m^2 to km^2
HR[,"est"] <- "km^2" %#% HR[, "est"] # converts m^2 to km^2
HR[,"high"] <- "km^2" %#% HR[,"high"] # converts m^2 to km^2
colnames(HR) <- c("HR_low", "HR_est", "HR_high")
# tau velocity
TAUVELOCITY <- ctmm::meta(FO_FITS, variable = "tau velocity", level = 0.95, units= FALSE)
TAUVELOCITY <- as.data.frame(TAUVELOCITY)
TAUVELOCITY <- TAUVELOCITY[-c(2,3),] # removes unnecessary rows
TAUVELOCITY[,"low"] <- "minutes" %#% TAUVELOCITY[,"low"] # converts seconds to minutes
TAUVELOCITY[,"est"] <- "minutes" %#% TAUVELOCITY[, "est"] # converts seconds to minutes
TAUVELOCITY[,"high"] <- "minutes" %#% TAUVELOCITY[,"high"] # converts seconds to minutes
colnames(TAUVELOCITY) <- c("τvelocity_low", "τvelocity_est", "τvelocity_high")
# combine all
META_total_FO <- cbind(SPEED, TAUVELOCITY, HR)
META_total_FO$Sex <- "Female"
META_total_FO$Status <- "Orphaned"


# male orphaned populations
# speed
SPEED <- ctmm::meta(MO_SPEED, variable = "speed", level = 0.95, units = FALSE)
SPEED <- as.data.frame(SPEED)
SPEED <- SPEED[-c(2,3),] # removes unnecessary rows
SPEED[,"low"] <- "km/day" %#% SPEED[,"low"] # converts m/s to km/day 
SPEED[,"est"] <- "km/day" %#% SPEED[, "est"] # converts m/s to km/day 
SPEED[,"high"] <- "km/day" %#% SPEED[,"high"] # converts m/s to km/day 
colnames(SPEED) <- c("speed_low", "speed_est", "speed_high")

# home range size 
HR <- ctmm::meta(MO_AKDE, variable = "area", level = 0.95, units = FALSE)
HR <- as.data.frame(HR)
HR <- HR[-c(2,3),] # removes unnecessary rows
HR[,"low"] <- "km^2" %#% HR[,"low"] # converts m^2 to km^2
HR[,"est"] <- "km^2" %#% HR[, "est"] # converts m^2 to km^2
HR[,"high"] <- "km^2" %#% HR[,"high"] # converts m^2 to km^2
colnames(HR) <- c("HR_low", "HR_est", "HR_high")
# tau velocity
TAUVELOCITY <- ctmm::meta(MO_FITS, variable = "tau velocity", level = 0.95, units= FALSE)
TAUVELOCITY <- as.data.frame(TAUVELOCITY)
TAUVELOCITY <- TAUVELOCITY[-c(2,3),] # removes unnecessary rows
TAUVELOCITY[,"low"] <- "minutes" %#% TAUVELOCITY[,"low"] # converts seconds to minutes
TAUVELOCITY[,"est"] <- "minutes" %#% TAUVELOCITY[, "est"] # converts seconds to minutes
TAUVELOCITY[,"high"] <- "minutes" %#% TAUVELOCITY[,"high"] # converts seconds to minutes
colnames(TAUVELOCITY) <- c("τvelocity_low", "τvelocity_est", "τvelocity_high")
# combine all
META_total_MO <- cbind(SPEED, TAUVELOCITY, HR)
META_total_MO$Sex <- "Male"
META_total_MO$Status <- "Orphaned"


# bind them all
META_df <- rbind(META_total_FO, META_total_MO)
# reorder IDs in alpbetical order
Movement_df <- Movement_df %>%
  arrange(ID)

# ensure IDs are characters for plotting
Movement_df$ID <- as.character(Movement_df$ID) 

# add sex information
Movement_df[Movement_df$ID %in% c("Bumpus", "Cate", "Elaine", "Makao", "Puji", "Annie", "Hannah", "Jane",
                                  "Margaret", "Maria", "Sheron", "Arya", "Dumbo_1", "Dumbo_2", "Colete", 
                                  "Heather", "Juju_2", "Mulan", "Rita", "Bella", "Nancy"), "Sex"] <- "Female"

Movement_df[Movement_df$ID %in% c("Anthony", "Christoffer", "Jackson", "Kyle", "Little_Rick", "Alexander",
                                  "Beto", "Larry", "Luigi", "Reid", "Rodolfo", "Thomas","Peter", "George", 
                                  "Dom", "Tim_3"), "Sex"] <- "Male"

# compare home range size between sexes
HR <- ggplot() +
  #total wild-raised
  geom_pointrange(data = subset(Movement_df, Status == "Orphaned" & Sex == "Female"), 
                  aes(x = HR_est, y = ID, xmin = HR_low, xmax = HR_high, color = "Female"), 
                  shape = 17) +
  geom_vline(xintercept = META_df[1,"HR_est"],
             color = "#E823A1", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[1,], aes(xmin = HR_low, xmax = HR_high, ymin = -Inf, ymax = Inf, color = "Female"), 
            color = NA, fill = "#E823A1", alpha = 0.15) +
  #total orphaned
  geom_pointrange(data = subset(Movement_df, Status == "Orphaned" & Sex == "Male"), 
                  aes(x = HR_est, y = ID, xmin = HR_low, xmax = HR_high, color = "Male"), 
                  shape = 15) +
  geom_vline(xintercept = META_df[2,"HR_est"], color = "#2919EC", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[2,], 
            aes(xmin = HR_low, xmax = HR_high, ymin = -Inf, ymax = Inf, color = "Male"), 
            color = NA, fill = "#2919EC", alpha = 0.1) +
  scale_color_manual("Sex", 
                     values = c("Female" = "#E823A1", "Male" = "#2919EC")) +
  labs(title = "a.", 
       y = "", 
       x = "log(Home Range (km^2))") +
  scale_x_continuous(trans = "log10") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "darkgray"), 
        plot.title = element_text(size = 12, family = "sans", hjust = 0.025),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 12),
        legend.text = element_text(size=12, family = "sans", face = "bold"),)  
#plot(HR)



# compare directional persistence between sexes
tauvelocity <- ggplot() +
  #total wild-raised
  geom_pointrange(data = subset(Movement_df, Status == "Orphaned" & Sex == "Female"), 
                  aes(x = τvelocity_est, y = ID, xmin = τvelocity_low, xmax = τvelocity_high, color = "Female"), 
                  shape = 17) +
  geom_vline(xintercept = META_df[1,"τvelocity_est"], 
             color = "#E823A1", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[1,], 
            aes(xmin = τvelocity_low, xmax = τvelocity_high, ymin = -Inf, ymax = Inf, color = "Female"), 
            color = NA, fill = "#E823A1", alpha = 0.15) +
  #total orphaned
  geom_pointrange(data = subset(Movement_df, Status == "Orphaned" & Sex == "Male"), 
                  aes(x = τvelocity_est, y = ID, xmin = τvelocity_low, xmax = τvelocity_high, color = "Male"), 
                  shape = 15) +
  geom_vline(xintercept = META_df[2,"τvelocity_est"], 
             color = "#2919EC", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[2,], 
            aes(xmin = τvelocity_low, xmax = τvelocity_high, ymin = -Inf, ymax = Inf, color = "Male"), 
            color = NA, fill = "#2919EC", alpha = 0.15) +
  scale_color_manual("Sex", values = c("Female" = "#E823A1", "Male" = "#2919EC")) +
  labs(title = "b.",
       y = "",
       x = "log(τvelocity (minutes))") +
  scale_x_continuous(trans = "log10") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "darkgray"), 
        plot.title = element_text(size = 12, family = "sans", hjust = 0.025), 
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 12),
        legend.text = element_text(size=12, family = "sans", face = "bold"),)   
#plot(tauvelocity)


# compare speed between sexes
speed <- ggplot() +
  #total wild-raised
  geom_pointrange(data = subset(Movement_df, Status == "Orphaned" & Sex == "Female" & is.finite(Mean_Speed_est) & is.finite(Mean_Speed_low) & is.finite(Mean_Speed_high)), 
                  aes(x = Mean_Speed_est, y = ID, xmin = Mean_Speed_low, xmax = Mean_Speed_high, color = "Female"), 
                  shape = 17) +
  geom_vline(xintercept = META_df[1,"speed_est"], 
             color = "#E823A1", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[1,], 
            aes(xmin = speed_low, xmax = speed_high, ymin = -Inf, ymax = Inf, color = "Female"), 
            color = NA, fill = "#E823A1", alpha = 0.15) +
  #total orphaned
  geom_pointrange(data = subset(Movement_df, Status == "Orphaned" & Sex == "Male" & is.finite(Mean_Speed_est) & is.finite(Mean_Speed_low) & is.finite(Mean_Speed_high)), 
                  aes(x = Mean_Speed_est, y = ID, xmin = Mean_Speed_low, xmax = Mean_Speed_high, color = "Male"), 
                  shape = 15) +
  geom_vline(xintercept = META_df[2,"speed_est"], 
             color = "#2919EC", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[2,], 
            aes(xmin = speed_low, xmax = speed_high, ymin = -Inf, ymax = Inf, color = "Male"), 
            color = NA, fill = "#2919EC", alpha = 0.15) +
  scale_color_manual("Sex", values = c("Female" = "#E823A1", "Male" = "#2919EC")) +
  labs(title = "c.", 
       y = "", 
       x = "Speed (km/day)") +
  #scale_x_continuous(trans = "log10") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "darkgray"), 
        plot.title = element_text(size = 12, family = "sans", hjust = 0.025),      
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 12),
        legend.text = element_text(size=12, family = "sans", face = "bold"),)   
#plot(speed)

# combine plots into one
movement_plot <- ggarrange(HR, tauvelocity, speed,
                           ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
# plot
plot(movement_plot)

# save final plot
ggsave("./FIGURES/Revised/Supplementary_Figure_3.png", width = 8.5, height = 5.5, units = "in", 
       dpi = 300, bg = "transparent")



