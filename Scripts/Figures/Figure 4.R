# load packages
library(ctmm) #movement
library(dplyr)
library(ggplot2) #making plots
library(ggpubr) #arranging plots in multiples
library(tidyr) #pivot_longer()
library(patchwork) # plotting final combo


# load results
load("./RESULTS/Window/Total_meta_RR_df.rda")

# plot ratios for home range size over time
HR_win <- ggplot() +
  geom_line(data = META_win_orphan, 
            aes(x = time_since_release, y = HR_est, color = "Ratio of Orphaned to Wild-Raised Population"), 
            size = 1) +
  geom_ribbon(data = META_win_orphan, 
              aes (x = time_since_release, ymin = HR_low, ymax = HR_high, fill = "Ratio of Orphaned to Wild-Raised Population"), 
              alpha = 0.1) +
  geom_hline(yintercept = 1, 
             col = "grey", 
             linetype = "dashed",
             size = 0.75) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c("Ratio of Orphaned to Wild-Raised Population" = "#23C3A8")) + 
  scale_color_manual(values = c("Ratio of Orphaned to Wild-Raised Population" = "#23C3A8")) + 
  labs(title = "a.", 
       y = "log(home range[orphan] / home range[wild-raised])", 
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



# plot ratios of speed over time
speed_win <- ggplot() +
  geom_line(data = META_win_orphan, 
            aes(x = time_since_release, y = speed_est, color = "Ratio of Orphaned to Wild-Raised Population"), 
            size = 1) +
  geom_ribbon(data = META_win_orphan, 
              aes (x = time_since_release, ymin = speed_low, ymax = speed_high, fill = "Ratio of Orphaned to Wild-Raised Population"), 
              alpha = 0.1) +
  geom_hline(yintercept = 1, 
             col = "grey", 
             linetype = "dashed",
             size = 0.75) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c("Ratio of Orphaned to Wild-Raised Population" = "#23C3A8")) + 
  scale_color_manual(values = c("Ratio of Orphaned to Wild-Raised Population" = "#23C3A8")) + 
  labs(title = "b.", y = "log(speed[orphan] / speed[wild-raised])",
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


# plot ratios of directional persistence over time
tauv_win <- ggplot() +
  geom_line(data = META_win_orphan, 
            aes(x = time_since_release, y = tauvelocity_est, color = "Ratio of Orphaned to Wild-Raised Population"), 
            size = 1) +
  geom_ribbon(data = META_win_orphan, 
              aes (x = time_since_release, ymin = tauvelocity_low, ymax = tauvelocity_high, fill = "Ratio of Orphaned to Wild-Raised Population"), 
              alpha = 0.1) +
  geom_hline(yintercept = 1, 
             col = "grey", 
             linetype = "dashed",
             size = 0.75) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c("Ratio of Orphaned to Wild-Raised Population" = "#23C3A8")) + 
  scale_color_manual(values = c("Ratio of Orphaned to Wild-Raised Population" = "#23C3A8")) + 
  labs(title = "c.", 
       y = "log(tau velocity[orphan] / tauvelocity[wild-raised])", 
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


# create layout for plotting
layout_design <- "
  A
  B
  C
"

# put plots together
HR_win + speed_win + tauv_win + plot_layout(design = layout_design) &   
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.title.position = "plot"
  )

# save plot
ggsave("./FIGURES/Revised/Figure_4.png", width = 8.5, height = 9.5, units = "in", 
       dpi = 300, bg = "transparent")


