# load packages
library(ctmm) #movement
library(dplyr)
library(ggplot2) #making plots
library(ggpubr) #arranging plots in multiples
library(tidyr) #drop_na function

# load RSFs
load("./RESULTS/RSFs_FOREST/OLD_0.055_1.5/RSF_wild.rda")

# separate RSFs based on age classes
SUBADULT_RSF <- RSF_wild[c("Anthony", "Kyle", "Little_Rick", "Reid")]

ADULT_RSF <- RSF_wild[c("Bumpus", "Cate", "Christoffer", "Elaine", "Jackson", "Makao", 
                         "Puji", "Alexander", "Annie", "Beto","Hannah", "Jane", "Larry", 
                         "Luigi", "Margaret", "Maria","Rodolfo", "Sheron", "Thomas")]

# get selection averages for each age class
ADULT_RSF_meta <- mean(ADULT_RSF)
SUBADULT_RSF_meta <- mean(SUBADULT_RSF)


# adult results ----
# create a dataframe to hold an individual's results
ADULT_RSF_df <- data.frame()

# extract summary output for an individual
df <- data.frame(t(summary(ADULT_RSF_meta, units = FALSE)[["CI"]]))
# if value is available in summary output, extract it into the individual's dataframe. 
# if value is not available in summary output, include an NA in the column
for(i in 1:1) {
  # Agriculture
  if ("Agriculture..1.Agriculture." %in% colnames(df)) {
    ADULT_RSF_df [i, "Agriculture_low"] <- df["low", "Agriculture..1.Agriculture."]
  }else{
    ADULT_RSF_df[i, "Agriculture_low"] <- NA
  }
  if ("Agriculture..1.Agriculture." %in% colnames(df)) {
    ADULT_RSF_df [i, "Agriculture_est"] <- df["est", "Agriculture..1.Agriculture."]
  }else{
    ADULT_RSF_df[i, "Agriculture_est"] <- NA
  }
  if ("Agriculture..1.Agriculture." %in% colnames(df)) {
    ADULT_RSF_df [i, "Agriculture_high"] <- df["high", "Agriculture..1.Agriculture."]
  }else{
    ADULT_RSF_df[i, "Agriculture_high"] <- NA
  }
  #Development
  if ("Development..1.Development." %in% colnames(df)) {
    ADULT_RSF_df [i, "Development_low"] <- df["low", "Development..1.Development."]
  }else{
    ADULT_RSF_df[i, "Development_low"] <- NA
  }
  if ("Development..1.Development." %in% colnames(df)) {
    ADULT_RSF_df [i, "Development_est"] <- df["est", "Development..1.Development."]
  }else{
    ADULT_RSF_df[i, "Development_est"] <- NA
  }
  if ("Development..1.Development." %in% colnames(df)) {
    ADULT_RSF_df [i, "Development_high"] <- df["high", "Development..1.Development."]
  }else{
    ADULT_RSF_df[i, "Development_high"] <- NA
  }
  # Forest
  if ("Forest..1.Forest." %in% colnames(df)) {
    ADULT_RSF_df [i, "Forest_low"] <- df["low", "Forest..1.Forest."]
  }else{
    ADULT_RSF_df[i, "Forest_low"] <- NA
  }
  if ("Forest..1.Forest." %in% colnames(df)) {
    ADULT_RSF_df [i, "Forest_est"] <- df["est", "Forest..1.Forest."]
  }else{
    ADULT_RSF_df[i, "Forest_est"] <- NA
  }
  if ("Forest..1.Forest." %in% colnames(df)) {
    ADULT_RSF_df [i, "Forest_high"] <- df["high", "Forest..1.Forest."]
  }else{
    ADULT_RSF_df[i, "Forst_high"] <- NA
  }
  
  # Pasture
  if ("Pasture..1.Pasture." %in% colnames(df)) {
    ADULT_RSF_df [i, "Pasture_low"] <- df["low", "Pasture..1.Pasture."]
  }else{
    ADULT_RSF_df[i, "Pasture_low"] <- NA
  }
  if ("Pasture..1.Pasture." %in% colnames(df)) {
    ADULT_RSF_df [i, "Pasture_est"] <- df["est", "Pasture..1.Pasture."]
  }else{
    ADULT_RSF_df[i, "Pasture_est"] <- NA
  }
  if ("Pasture..1.Pasture." %in% colnames(df)) {
    ADULT_RSF_df [i, "Pasture_high"] <- df["high", "Pasture..1.Pasture."]
  }else{
    ADULT_RSF_df[i, "Pasture_high"] <- NA
  }
  # Water
  if ("Water..1.Water." %in% colnames(df)) {
    ADULT_RSF_df [i, "Water_low"] <- df["low", "Water..1.Water."]
  }else{
    ADULT_RSF_df[i, "Water_low"] <- NA
  }
  if ("Water..1.Water." %in% colnames(df)) {
    ADULT_RSF_df [i, "Water_est"] <- df["est", "Water..1.Water."]
  }else{
    ADULT_RSF_df[i, "Water_est"] <- NA
  }
  if ("Water..1.Water." %in% colnames(df)) {
    ADULT_RSF_df [i, "Water_high"] <- df["high", "Water..1.Water."]
  }else{
    ADULT_RSF_df[i, "Water_high"] <- NA
  }} # end of inner loop




# subadult results ----

# create a dataframe to hold an individual's results
SUBADULT_RSF_df <- data.frame()

# extract summary output for an individual
df <- data.frame(t(summary(SUBADULT_RSF_meta, units = FALSE)[["CI"]]))
# if value is available in summary output, extract it into the individual's dataframe. 
# if value is not available in summary output, include an NA in the column
for(i in 1:1) {
  # Agriculture
  if ("Agriculture..1.Agriculture." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Agriculture_low"] <- df["low", "Agriculture..1.Agriculture."]
  }else{
    SUBADULT_RSF_df[i, "Agriculture_low"] <- NA
  }
  if ("Agriculture..1.Agriculture." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Agriculture_est"] <- df["est", "Agriculture..1.Agriculture."]
  }else{
    SUBADULT_RSF_df[i, "Agriculture_est"] <- NA
  }
  if ("Agriculture..1.Agriculture." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Agriculture_high"] <- df["high", "Agriculture..1.Agriculture."]
  }else{
    SUBADULT_RSF_df[i, "Agriculture_high"] <- NA
  }
  # Development
  if ("Development..1.Development." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Development_low"] <- df["low", "Development..1.Development."]
  }else{
    SUBADULT_RSF_df[i, "Development_low"] <- NA
  }
  if ("Development..1.Development." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Development_est"] <- df["est", "Development..1.Development."]
  }else{
    SUBADULT_RSF_df[i, "Development_est"] <- NA
  }
  if ("Development..1.Development." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Development_high"] <- df["high", "Development..1.Development."]
  }else{
    SUBADULT_RSF_df[i, "Development_high"] <- NA
  }
  # Forest
  if ("Forest..1.Forest." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Forest_low"] <- df["low", "Forest..1.Forest."]
  }else{
    SUBADULT_RSF_df[i, "Forest_low"] <- NA
  }
  if ("Forest..1.Forest." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Forest_est"] <- df["est", "Forest..1.Forest."]
  }else{
    SUBADULT_RSF_df[i, "Forest_est"] <- NA
  }
  if ("Forest..1.Forest." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Forest_high"] <- df["high", "Forest..1.Forest."]
  }else{
    SUBADULT_RSF_df[i, "Forst_high"] <- NA
  }
  
  # Pasture
  if ("Pasture..1.Pasture." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Pasture_low"] <- df["low", "Pasture..1.Pasture."]
  }else{
    SUBADULT_RSF_df[i, "Pasture_low"] <- NA
  }
  if ("Pasture..1.Pasture." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Pasture_est"] <- df["est", "Pasture..1.Pasture."]
  }else{
    SUBADULT_RSF_df[i, "Pasture_est"] <- NA
  }
  if ("Pasture..1.Pasture." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Pasture_high"] <- df["high", "Pasture..1.Pasture."]
  }else{
    SUBADULT_RSF_df[i, "Pasture_high"] <- NA
  }
  # Water
  if ("Water..1.Water." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Water_low"] <- df["low", "Water..1.Water."]
  }else{
    SUBADULT_RSF_df[i, "Water_low"] <- NA
  }
  if ("Water..1.Water." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Water_est"] <- df["est", "Water..1.Water."]
  }else{
    SUBADULT_RSF_df[i, "Water_est"] <- NA
  }
  if ("Water..1.Water." %in% colnames(df)) {
    SUBADULT_RSF_df [i, "Water_high"] <- df["high", "Water..1.Water."]
  }else{
    SUBADULT_RSF_df[i, "Water_high"] <- NA
  }} # end of inner loop


# add age class information to each dataframe
SUBADULT_RSF_df$Age <- "Subadult"
ADULT_RSF_df$Age <- "Adult"

# bind results together
mean_all <- rbind(SUBADULT_RSF_df, ADULT_RSF_df)  

# plot selection results
mean_all %>% 
  select(Age, contains("_low")) %>% 
  pivot_longer(cols = -c(Age), names_to = "type", values_to = "low") %>% 
  mutate(type = gsub("_low", "", type)) %>% 
  full_join(mean_all %>% 
              select(Age, contains("_est")) %>% 
              pivot_longer(cols = -c(Age), names_to = "type", values_to = "est") %>% 
              mutate(type = gsub("_est", "", type)), 
            by = c("Age", "type")) %>% 
  full_join(mean_all %>% 
              select(Age, contains("_high")) %>% 
              pivot_longer(cols = -c(Age), names_to = "type", values_to = "high") %>% 
              mutate(type = gsub("_high", "", type)), 
            by = c("Age", "type")) %>% 
  ggplot() +
  geom_errorbar(aes(x = est, xmin = low, xmax = high, y = type, color = Age), width = 0.5, linewidth = 2.5,  position = position_dodge(width = 0.5)) +
  scale_color_manual("Age", values = c("Subadult" = "#DA1C1C", "Adult" = "#8732D2")) +
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



ggsave("./FIGURES/Revised/Figure_S2.png", width = 8.5, height = 5.5, units = "in", 
       dpi = 300, bg = "transparent")




