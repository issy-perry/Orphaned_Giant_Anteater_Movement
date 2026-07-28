# load packages
library(ctmm)
library(ggplot2)
library(ggpubr)
library(dplyr)

# load data
load("./RESULTS/META_average_RR_df.rda")
load("./RESULTS/Movement_model_RR_df.rda")

# reorder IDs in alpbetical order
Movement_df <- Movement_df %>%
  arrange(ID)
# convert ID values to numeric for plotting
Movement_df$ID <- as.character(Movement_df$ID) 

# home range comparison plot ----
HR <- ggplot() +
  #total wild-raised
  geom_pointrange(data = subset(Movement_df, Status == "Wild-raised"), 
                  aes(x = HR_est, y = ID, xmin = HR_low, xmax = HR_high, color = "Wild-raised"), 
                  shape = 17) +
  geom_vline(xintercept = META_df[2,"HR_est"],
             color = "#EA8109", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[2,], aes(xmin = HR_low, xmax = HR_high, ymin = -Inf, ymax = Inf, color = "Wild-raised"), 
            color = NA, fill = "#EA8109", alpha = 0.15) +
  #total orphaned
  geom_pointrange(data = subset(Movement_df, Status == "Orphaned"), 
                  aes(x = HR_est, y = ID, xmin = HR_low, xmax = HR_high, color = "Orphaned"), 
                  shape = 15) +
  geom_vline(xintercept = META_df[1,"HR_est"], color = "#23C3A8", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[1,], 
            aes(xmin = HR_low, xmax = HR_high, ymin = -Inf, ymax = Inf, color = "Orphaned"), 
            color = NA, fill = "#23C3A8", alpha = 0.1) +
  scale_color_manual("Individual Type", 
                     values = c("Orphaned" = "#23C3A8", "Wild-raised" = "#EA8109")) +
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

# tau velocity comparison plot ----
tauvelocity <- ggplot() +
  #total wild-raised
  geom_pointrange(data = subset(Movement_df, Status == "Wild-raised"), 
                  aes(x = τvelocity_est, y = ID, xmin = τvelocity_low, xmax = τvelocity_high, color = "Wild-raised"), 
                  shape = 17) +
  geom_vline(xintercept = META_df[2,"τvelocity_est"], 
             color = "#EA8109", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[2,], 
            aes(xmin = τvelocity_low, xmax = τvelocity_high, ymin = -Inf, ymax = Inf, color = "Wild-raised"), 
            color = NA, fill = "#EA8109", alpha = 0.15) +
  #total orphaned
  geom_pointrange(data = subset(Movement_df, Status == "Orphaned"), 
                  aes(x = τvelocity_est, y = ID, xmin = τvelocity_low, xmax = τvelocity_high, color = "Orphaned"), 
                  shape = 15) +
  geom_vline(xintercept = META_df[1,"τvelocity_est"], 
             color = "#23C3A8", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[1,], 
            aes(xmin = τvelocity_low, xmax = τvelocity_high, ymin = -Inf, ymax = Inf, color = "Orphaned"), 
            color = NA, fill = "#23C3A8", alpha = 0.15) +
  scale_color_manual("Individual Type", values = c("Orphaned" = "#23C3A8", "Wild-raised" = "#EA8109")) +
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


# speed comparison plot ----
speed <- ggplot() +
  #total wild-raised
  geom_pointrange(data = subset(Movement_df, Status == "Wild-raised" & is.finite(Mean_Speed_est) & is.finite(Mean_Speed_low) & is.finite(Mean_Speed_high)), 
                  aes(x = Mean_Speed_est, y = ID, xmin = Mean_Speed_low, xmax = Mean_Speed_high, color = "Wild-raised"), 
                  shape = 17) +
  geom_vline(xintercept = META_df[2,"speed_est"], 
             color = "#EA8109", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[2,], 
            aes(xmin = speed_low, xmax = speed_high, ymin = -Inf, ymax = Inf, color = "Wild-raised"), 
            color = NA, fill = "#EA8109", alpha = 0.15) +
  #total orphaned
  geom_pointrange(data = subset(Movement_df, Status == "Orphaned" & is.finite(Mean_Speed_est) & is.finite(Mean_Speed_low) & is.finite(Mean_Speed_high)), 
                  aes(x = Mean_Speed_est, y = ID, xmin = Mean_Speed_low, xmax = Mean_Speed_high, color = "Orphaned"), 
                  shape = 15) +
  geom_vline(xintercept = META_df[1,"speed_est"], 
             color = "#23C3A8", lty = "dashed", lwd = 1) + #area_est
  geom_rect(data = META_df[1,], 
            aes(xmin = speed_low, xmax = speed_high, ymin = -Inf, ymax = Inf, color = "Orphaned"), 
            color = NA, fill = "#23C3A8", alpha = 0.15) +
  scale_color_manual("Individual Type", values = c("Orphaned" = "#23C3A8", "Wild-raised" = "#EA8109")) +
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


movement_plot <- ggarrange(HR, tauvelocity, speed,
                           ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")

plot(movement_plot)

ggsave("./FIGURES/Revised/Figure_2.png", width = 8.5, height = 5.5, units = "in", 
       dpi = 300, bg = "transparent")
