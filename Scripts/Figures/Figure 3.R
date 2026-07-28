#load packages
library(ctmm) #movement
library(dplyr)
library(ggplot2) #making plots
library(ggpubr) #arranging plots in multiples
library(tidyr) #drop_na function

load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/Mean_RSF_all_df.rda")
mean_all$Residence <- NA
O_W <- mean_all %>% 
  select(ID, Residence, contains("_low")) %>% 
  pivot_longer(cols = -c(ID, Residence), names_to = "type", values_to = "low") %>% 
  mutate(type = gsub("_low", "", type)) %>% 
  full_join(mean_all %>% 
              select(ID, Residence, contains("_est")) %>% 
              pivot_longer(cols = -c(ID, Residence), names_to = "type", values_to = "est") %>% 
              mutate(type = gsub("_est", "", type)), 
            by = c("ID", "Residence", "type")) %>% 
  full_join(mean_all %>% 
              select(ID, Residence, contains("_high")) %>% 
              pivot_longer(cols = -c(ID, Residence), names_to = "type", values_to = "high") %>% 
              mutate(type = gsub("_high", "", type)), 
            by = c("ID", "Residence", "type")) %>% 
  
  #23C3A8
  
  ggplot() +
  geom_errorbar(aes(x = est, xmin = low, xmax = high, y = type, color = ID), width = 0.5, linewidth = 2.5,  position = position_dodge(width = 0.5)) +
  scale_color_manual("Population", values = c("Orphaned" = "#23C3A8", "Wild-raised" = "#EA8109")) +
  geom_vline(xintercept = 0, col = "grey70", linetype = "dashed", linewidth = 1) +
  labs(y = "", 
       x = "Poisson Regression Coefficients", 
       title = "") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "darkgray"), 
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))

plot(O_W)




ggsave("./FIGURES/Revised/Figure_3_with_development.png", width = 8.5, height = 5.5, units = "in", 
       dpi = 300, bg = "transparent")








# without development ------

mean_all[,c("Development_low", "Development_est", "Development_high")] <- NULL
O_W <- mean_all %>% 
  select(ID, Residence, contains("_low")) %>% 
  pivot_longer(cols = -c(ID, Residence), names_to = "type", values_to = "low") %>% 
  mutate(type = gsub("_low", "", type)) %>% 
  full_join(mean_all %>% 
              select(ID, Residence, contains("_est")) %>% 
              pivot_longer(cols = -c(ID, Residence), names_to = "type", values_to = "est") %>% 
              mutate(type = gsub("_est", "", type)), 
            by = c("ID", "Residence", "type")) %>% 
  full_join(mean_all %>% 
              select(ID, Residence, contains("_high")) %>% 
              pivot_longer(cols = -c(ID, Residence), names_to = "type", values_to = "high") %>% 
              mutate(type = gsub("_high", "", type)), 
            by = c("ID", "Residence", "type")) %>% 
  
  
  
  ggplot() +
  geom_errorbar(aes(x = est, xmin = low, xmax = high, y = type, color = ID), width = 0.5, linewidth = 2.5,  position = position_dodge(width = 0.5)) +
  scale_color_manual("Population", values = c("Orphaned" = "#23C3A8", "Wild-raised" = "#EA8109")) +
  geom_vline(xintercept = 0, col = "grey70", linetype = "dashed", linewidth = 1) +
  labs(y = "", 
       x = "Poisson Regression Coefficients",
       title = "") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "darkgray"), 
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size=12, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))

plot(O_W)




ggsave("./FIGURES/Revised/Figure_3_no_development.png", width = 8.5, height = 5.5, units = "in", 
       dpi = 300, bg = "transparent")

