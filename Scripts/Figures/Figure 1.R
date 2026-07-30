# load packages
library(ggplot2)
library(ctmm)
library(terra)
library(sf)
library(tidyterra)
library(gridExtra)
library(rphylopic)
library(ggspatial)
library(tidyr) #pivot_longer()
library(cowplot) # plotting 
library(patchwork)

# load raster data
# not using pre-written script since we have different extents for plotting
land_cover_2024 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2024.tif")
land_cover_2023 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2023.tif")
land_cover_2022 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2022.tif")
land_cover_2021 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2021.tif")
land_cover_2020 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2020.tif")
land_cover_2019 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2019.tif")
land_cover_2018 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2018.tif")
land_cover_2017 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2017.tif")


# load tel data
load("./DATA/Wild_raised_tel_data/Data_telemetry_RR.rda") # wild-raised 
load("./DATA/Orphaned_tel_data/Data_telemetry_RR.rda") # orphaned 

# telemetry overlayed on raster ------------------------------------------------------------------------------------------------------------
# separate orphaned population into different study sites for better plotting
DATA_orphan_b <- DATA_orphan[c("Arya", "Dumbo_1", "Dumbo_2", "Tim_2", "Colete", "Heather", 
                               "Mulan", "Tim_3", "Dom", "Rita", "Bella", "George")]

DATA_orphan_c <- DATA_orphan[c("Juju_2", "Peter", "Nancy")]

# create new extents
ew <- ext(-53.9, -53.45, -21.79, -20.98) 
eb <- ext(-48.34, -47.72, -19.23, -18.7)
ec <- ext(-48.81, -48.69, -18.45, -18.34)

# crop rasters
cover_2017 <- crop(land_cover_2017, ew)
cover_2018 <- crop(land_cover_2018, ew)
cover_2019 <- crop(land_cover_2019, eb)
cover_2019c <- crop(land_cover_2019, ec)
cover_2020 <- crop(land_cover_2020, eb)
cover_2020c <- crop(land_cover_2020, ec)
cover_2021 <- crop(land_cover_2021, eb)
cover_2021c <- crop(land_cover_2021, ec)
cover_2022 <- crop(land_cover_2022, eb)
cover_2022c <- crop(land_cover_2022, ec)
cover_2023 <- crop(land_cover_2023, eb)
cover_2023c <- crop(land_cover_2023, ec)
cover_2024 <- crop(land_cover_2024, eb)
cover_2024c <- crop(land_cover_2024, ec)
# clear environment space
rm(land_cover_2017, land_cover_2018, land_cover_2019, land_cover_2020, land_cover_2021, land_cover_2022, land_cover_2023, land_cover_2024)
gc()

# plot study site a
land_types <- cover_2018
DATA_wild_df <- do.call(rbind.data.frame, DATA_wild)
DATA_wild_sf <- st_as_sf(DATA_wild_df, coords = c("longitude", "latitude"), 
                         crs = crs(land_types))
DATA_wild_ext <- ext(DATA_wild_sf)
#SpatExtent : -53.79925, -53.474271, -21.772054, -21.08363 (xmin, xmax, ymin, ymax)
total_ext <- ext(-53.9, -53.45, -21.79, -20.98) #create extent barriers so that we can crop a certain amount out for orphaned anteaters

land_cropped <- crop(land_types,
                     total_ext,
                     snap = "out")

land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)


wild_plot <- ggplot() +
  geom_spatraster(data = land_cropped, maxcell = 5e+07,
                  aes(fill = brazil_coverage_2018)) +
  geom_sf(data = DATA_wild_sf, size = 1.45, alpha = 0.6, col = "black", shape = 16) +
  scale_fill_manual(breaks = c(94,95,96,97,98),
                    labels = c("Agriculture","Development","Forest","Pasture", "Water"),
                    values = c("#B7E4B9", "#928680", "#1A583F", "#F9F5C1", "#94BEE0"),
                    name = "Land Class") +
  #theme_bw() +
  labs(title = "a.") +
  theme(plot.title = element_text(size = 12, hjust = 0.35),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        #legend.position.inside = c(0.28,0.12),
        legend.position = "none",
        legend.title = element_text(size=14, family = "sans", face = "bold", hjust = 0.5),
        legend.text = element_text(size=12, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        #strip.background=element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggspatial::annotation_scale(
    location = "tr",
    pad_y = unit(-0.007, "cm"),
    bar_cols = c("white", "#404040"),
    text_family = "sans") +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.3, "cm"), pad_y = unit(0.4, "cm"),
    height = unit(0.5, "cm"), width = unit(0.5, "cm"),
    style = north_arrow_orienteering(
      text_face = "bold",
      text_size = 10,
      text_col = "black",
      fill = c("white", "#404040"))) 




# plot study site b
land_types <- cover_2021
DATA_orphan_df <- do.call(rbind.data.frame, DATA_orphan_b)
DATA_orphan_b_sf <- st_as_sf(DATA_orphan_df, coords = c("longitude", "latitude"), 
                           crs = crs(land_types))
DATA_orphan_ext <- ext(DATA_orphan_b_sf)
#SpatExtent :  -48.319349, -47.742377, -19.210794, -18.725847 (xmin, xmax, ymin, ymax)
#separate into multiple study sites so that we can actually see the telemetry data
total_ext <- ext(-48.34, -47.72, -19.23, -18.7)

land_cropped <- crop(land_types,
                     total_ext,
                     snap = "out")

#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)


orphan_b_plot <- ggplot() +
  geom_spatraster(data = land_cropped, maxcell = 5e+07,
                  aes(fill = brazil_coverage_2021)) +
  geom_sf(data = DATA_orphan_b_sf, size = 1.5, alpha = 0.6, col = "black", shape = 16) +
  scale_fill_manual(breaks = c(94,95,96,97,98),
                    labels = c("Agriculture","Development","Forest","Pasture", "Water"),
                    values = c("#B7E4B9", "#928680", "#1A583F", "#F9F5C1", "#94BEE0"),
                    name = "Land Class") +
  labs(title = "b.") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.05),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        #legend.position.inside = c(0.28,0.12),
        legend.position = "none",
        legend.title = element_text(size=12, family = "sans", face = "bold", hjust = 0.5), #element_text(size=12, family = "sans", face = "bold", hjust = 0.5
        legend.text = element_text(size=11, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        #strip.background=element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggspatial::annotation_scale(
    location = "tr",
    pad_y = unit(-0.007, "cm"),
    bar_cols = c("white", "#404040"),
    text_family = "sans") +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.8, "cm"), pad_y = unit(0.5, "cm"),
    height = unit(0.5, "cm"), width = unit(0.5, "cm"),
    style = north_arrow_orienteering(
      text_face = "bold",
      text_size = 10,
      text_col = "black",
      fill = c("white", "#404040"))) 


# plot study site c
land_types <- cover_2021c
DATA_orphan_df <- do.call(rbind.data.frame, DATA_orphan_c)
DATA_orphan_c_sf <- st_as_sf(DATA_orphan_df, coords = c("longitude", "latitude"), 
                           crs = crs(land_types))
DATA_orphan_ext <- ext(DATA_orphan_c_sf)
#SpatExtent : -48.79525, -48.714981, -18.432612, -18.362662 (xmin, xmax, ymin, ymax)
#separate into multiple study sites so that we can actually see the telemetry data
total_ext <- ext(-48.81, -48.69, -18.45, -18.34)

land_cropped <- crop(land_types,
                     total_ext,
                     snap = "out")

#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)


orphan_c_plot <- ggplot() +
  geom_spatraster(data = land_cropped, maxcell = 5e+07,
                  aes(fill = brazil_coverage_2021)) +
  geom_sf(data = DATA_orphan_c_sf, size = 1.5, alpha = 0.6, col = "black", shape = 16) +
  scale_fill_manual(breaks = c(94,95,96,97,98),
                    labels = c("Agriculture","Development","Forest","Pasture", "Water"),
                    values = c("#B7E4B9", "#928680", "#1A583F", "#F9F5C1", "#94BEE0"),
                    name = "Land Class") +
  labs(title = "c.") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.035),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        #legend.position.inside = c(0.28,0.12),
        legend.position = "right",
        legend.title = element_text(size=12, family = "sans", face = "bold", hjust = 0.5), #element_text(size=12, family = "sans", face = "bold", hjust = 0.5
        legend.text = element_text(size=11, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        #strip.background=element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggspatial::annotation_scale(
    location = "tr",
    pad_y = unit(-0.007, "cm"),
    bar_cols = c("white", "#404040"),
    text_family = "sans") +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.8, "cm"), pad_y = unit(0.5, "cm"),
    height = unit(0.5, "cm"), width = unit(0.5, "cm"),
    style = north_arrow_orienteering(
      text_face = "bold",
      text_size = 10,
      text_col = "black",
      fill = c("white", "#404040"))) 

#plot(orphan_plot)





# proportions of land coverage -------------------------------------------------------------------------------
# wild-rasied sites 
land_types <- cover_2017
DATA_wild_df <- do.call(rbind.data.frame, DATA_wild)
DATA_wild_sf <- st_as_sf(DATA_wild_df, coords = c("longitude", "latitude"), 
                         crs = crs(land_types))
DATA_wild_ext <- ext(DATA_wild_sf)
land_cropped <- crop(land_types,
                     DATA_wild_ext,
                     snap = "out")

#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_17 <- data.frame()
prop_17[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_17[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_17[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_17[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_17[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

#2018
land_types <- cover_2018
land_cropped <- crop(land_types,
                     DATA_wild_ext,
                     snap = "out")

#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_18 <- data.frame()
prop_18[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_18[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_18[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_18[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_18[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

# add information about years
prop_17$Year <- "2017"
prop_18$Year <- "2018"
# bind dataframes together
prop_w <- rbind(prop_17, prop_18)
prop_w$Site <- "Wild-raised"



# orphaned sites
#2019 
land_cropped <- cover_2019
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_19 <- data.frame()
prop_19[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_19[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_19[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_19[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_19[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

#2020 
land_cropped <- cover_2020
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_20 <- data.frame()
prop_20[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_20[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_20[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_20[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_20[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster



#2021 
land_cropped <- cover_2021
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_21 <- data.frame()
prop_21[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_21[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_21[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_21[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_21[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

#2022 
land_cropped <- cover_2022
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_22 <- data.frame()
prop_22[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_22[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_22[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_22[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_22[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

#2023 
land_cropped <- cover_2023
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_23 <- data.frame()
prop_23[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_23[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_23[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_23[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_23[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

#2024 
land_cropped <- cover_2024
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_24 <- data.frame()
prop_24[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_24[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_24[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_24[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_24[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

# add information about years
prop_19$Year <- "2019"
prop_20$Year <- "2020"
prop_21$Year <- "2021"
prop_22$Year <- "2022"
prop_23$Year <- "2023"
prop_24$Year <- "2024"




# 2019 c.
land_cropped <- cover_2019c
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_19c <- data.frame()
prop_19c[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_19c[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_19c[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_19c[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_19c[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

#2020 c.
land_cropped <- cover_2020c
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_20c <- data.frame()
prop_20c[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_20c[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_20c[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_20c[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_20c[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

#2021 c.
land_cropped <- cover_2021c
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_21c <- data.frame()
prop_21c[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_21c[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_21c[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_21c[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_21c[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

#2022 c.
land_cropped <- cover_2022c
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_22c <- data.frame()
prop_22c[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_22c[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_22c[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_22c[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_22c[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

#2023 c.
land_cropped <- cover_2023c
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_23c <- data.frame()
prop_23c[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_23c[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_23c[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_23c[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_23c[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

#2024 c.
land_cropped <- cover_2024c
#Process the land class raster for easy plotting
land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29,9), 96, land_cropped) #97 = Native_forest
land_cropped <- ifel(land_cropped %in% c(12,15), 97, land_cropped) #98 = Pasture
land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 94, land_cropped) # 94 = Agriculture
land_cropped <- ifel(land_cropped %in% c(24,25,50,75,30), 95, land_cropped) # 95 = Development
land_cropped <- ifel(land_cropped %in% c(11,26,33,31), 98, land_cropped) # 99 = water
land_cropped <- ifel(land_cropped %in% 0:93, 93, land_cropped) # 93 = other
land_cropped <- as.factor(land_cropped)

prop_24c <- data.frame()
prop_24c[1, "Agriculture"] <- prop.table(table(as.vector(land_cropped))) ["94"] #proportion of agriculture land type in raster
prop_24c[1, "Development"] <- prop.table(table(as.vector(land_cropped))) ["95"] #proportion of development land type in raster
prop_24c[1, "Forest"] <- prop.table(table(as.vector(land_cropped))) ["96"] #proportion of native forest land type in raster
prop_24c[1, "Pasture"] <- prop.table(table(as.vector(land_cropped))) ["97"] #proportion of pasture land type in raster
prop_24c[1, "Water"] <- prop.table(table(as.vector(land_cropped))) ["98"] #proportion of agriculture land type in raster

# add information about years
prop_19c$Year <- "2019"
prop_20c$Year <- "2020"
prop_21c$Year <- "2021"
prop_22c$Year <- "2022"
prop_23c$Year <- "2023"
prop_24c$Year <- "2024"
# bind dataframes together
prop_o <- rbind(prop_19, prop_20, prop_21, prop_22, prop_23, prop_24,
                prop_19c, prop_20c, prop_21c, prop_22c, prop_23c, prop_24c)
prop_o$Site <- "Orphaned"

#means of proportions across years wild-raised site
mean_w <- data.frame()
mean_w[1, "Agriculture"] <- mean(prop_w[["Agriculture"]])
mean_w[1, "Development"] <- mean(prop_w[["Development"]])
mean_w[1, "Forest"] <- mean(prop_w[["Forest"]])
mean_w[1, "Pasture"] <- mean(prop_w[["Pasture"]])
mean_w[1, "Water"] <- mean(prop_w[["Water"]])

mean_w$Year <- "Average_W"
mean_w$Site <- "Wild-raised"

#means of proportions across years orphan site
mean_o <- data.frame()
mean_o[1, "Agriculture"] <- mean(prop_o[["Agriculture"]])
mean_o[1, "Development"] <- mean(prop_o[["Development"]])
mean_o[1, "Forest"] <- mean(prop_o[["Forest"]])
mean_o[1, "Pasture"] <- mean(prop_o[["Pasture"]])
mean_o[1, "Water"] <- mean(prop_o[["Water"]])

mean_o$Year <- "Average_O"
mean_o$Site <- "Orphaned"

#combine dataframes
prop_w <- rbind(prop_w, mean_w)
prop_o <- rbind(prop_o, mean_o)
prop_land <-rbind(prop_w, prop_o)

# convert proportions to percentages
prop_land[, 1:5] <- prop_land[, 1:5] * 100

# free environment space
rm(prop_17, prop_18, mean_w, prop_19, prop_20, prop_21, prop_22, prop_23, prop_24, mean_o)
save(prop_land, file = "./RESULTS/Land_Use/Total_proportion_df_RR.rda")

#pivot to make plotting easier
prop_land_piv <- prop_land %>% pivot_longer(cols=c("Water", "Agriculture", "Development", "Pasture", "Forest"),
                                            names_to = "Land_Type",
                                            values_to = "Proportion")
save(prop_land_piv, file = "./RESULTS/Land_Use/Total_proportion_pivoted_RR.rda")
load("./RESULTS/Land_Use/Total_proportion_pivoted_RR.rda")

#remove years for now just to get average
prop_land_piv <- subset(prop_land_piv, Year != 2017)
prop_land_piv <- subset(prop_land_piv, Year != 2018)
prop_land_piv <- subset(prop_land_piv, Year != 2019)
prop_land_piv <- subset(prop_land_piv, Year != 2020)
prop_land_piv <- subset(prop_land_piv, Year != 2021)
prop_land_piv <- subset(prop_land_piv, Year != 2022)
prop_land_piv <- subset(prop_land_piv, Year != 2023)
prop_land_piv <- subset(prop_land_piv, Year != 2024)

#plot proportions
prop_land_piv$Site <- factor(prop_land_piv$Site, levels = c("Wild-raised", "Orphaned"))


land_plot <- ggplot(prop_land_piv, aes(x = Land_Type, y = Proportion, fill = Site)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual("Population", 
                    values = c("Wild-raised" = "#F5AD5B", "Orphaned" = "#79DECD"),
                    name = "Population") +
  labs(title = "d.", x = "", y = "Percentage of Coverage") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.065),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA), 
        legend.position = "right", 
        legend.title = element_text(size=12, family = "sans", face = "bold", hjust = 0.5),
        legend.text = element_text(size = 11, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        #strip.text = element_text(face = "bold", size = 12, family = "sans"),
        axis.line = element_line(color = "darkgray", linewidth = 0.5),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 12)) 

layout_design <- "
  ABC
  DDD
"

# 3. Assemble them together
wild_plot + orphan_b_plot + orphan_c_plot + land_plot + plot_layout(design = layout_design) &   
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.title.position = "plot"
  )

ggsave("./FIGURES/Revised/Figure_1_RR.png", width = 8.5, height = 6, units = "in", 
       dpi = 300, bg = "transparent")


