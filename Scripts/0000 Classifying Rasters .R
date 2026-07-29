# 0000 Classifying Rasters
# This script includes the following:
#     1. Call script for downloading raster datasets from Mapbiomas
#     2. Create new categorical rasters for each habitat (i.e. Agriculture Spatraster has 0s meaning no agriculture and 1s meaning agriculture is present)
#     3. Combine rasters into one list

# load packages
library(terra)
library(raster)

# create function for establishing land classes
`%notin%` <- Negate(`%in%`)

# call raster downloading script to run
source("./SCRIPTS/000 Downloading Land Coverage Rasters.R")

# 2024 ----
Agriculture <- cover_2024
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- FALSE
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- TRUE
Agriculture <- raster(Agriculture)

Development <- cover_2024
Development[Development %notin% c(24,25,30)] <- FALSE
Development[Development %in% c(24,25,30)] <- TRUE
Development <- raster(Development)

Forest <- cover_2024
Forest[Forest %notin% c(1,3,4,5,6,49,29,9)] <- FALSE
Forest[Forest %in% c(1,3,4,5,6,49,29,9)] <- TRUE
Forest <- raster(Forest)

Pasture <- cover_2024
Pasture[Pasture %notin% c(12,15)] <- FALSE
Pasture[Pasture %in% c(12,15)] <- TRUE
Pasture <- raster(Pasture)

Water <- cover_2024
Water[Water %notin% c(11,26,33)] <- FALSE
Water[Water %in% c(11,26,33)] <- TRUE
Water <- raster(Water)

covers_2024 <- list(Agriculture = Agriculture,
               Development = Development,
               Forest = Forest,
               Pasture = Pasture,
               Water = Water)

# free environment space
rm(Native_Forest, Forestry, Pasture, Agriculture, Development, Water, cover_2024)
gc()


# 2023 ----
Agriculture <- cover_2023
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- FALSE
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- TRUE
Agriculture <- raster(Agriculture)

Development <- cover_2023
Development[Development %notin% c(24,25,30)] <- FALSE
Development[Development %in% c(24,25,30)] <- TRUE
Development <- raster(Development)

Forest <- cover_2023
Forest[Forest %notin% c(1,3,4,5,6,49,29,9)] <- FALSE
Forest[Forest %in% c(1,3,4,5,6,49,29,9)] <- TRUE
Forest <- raster(Forest)

Pasture <- cover_2023
Pasture[Pasture %notin% c(12,15)] <- FALSE
Pasture[Pasture %in% c(12,15)] <- TRUE
Pasture <- raster(Pasture)

Water <- cover_2023
Water[Water %notin% c(11,26,33)] <- FALSE
Water[Water %in% c(11,26,33)] <- TRUE
Water <- raster(Water)

covers_2023 <- list(Agriculture = Agriculture,
               Development = Development,
               Forest = Forest,
               Pasture = Pasture,
               Water = Water)

# free environment space
rm(Native_Forest, Forestry, Pasture, Agriculture, Development, Water, cover_2023)
gc()


# 2022 ----
Agriculture <- cover_2022
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- FALSE
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- TRUE
Agriculture <- raster(Agriculture)

Development <- cover_2022
Development[Development %notin% c(24,25,30)] <- FALSE
Development[Development %in% c(24,25,30)] <- TRUE
Development <- raster(Development)

Forest <- cover_2022
Forest[Forest %notin% c(1,3,4,5,6,49,29,9)] <- FALSE
Forest[Forest %in% c(1,3,4,5,6,49,29,9)] <- TRUE
Forest <- raster(Forest)

Pasture <- cover_2022
Pasture[Pasture %notin% c(12,15)] <- FALSE
Pasture[Pasture %in% c(12,15)] <- TRUE
Pasture <- raster(Pasture)

Water <- cover_2022
Water[Water %notin% c(11,26,33)] <- FALSE
Water[Water %in% c(11,26,33)] <- TRUE
Water <- raster(Water)

covers_2022 <- list(Agriculture = Agriculture,
                    Development = Development,
                    Forest = Forest,
                    Pasture = Pasture,
                    Water = Water)

# free environment space
rm(Native_Forest, Forestry, Pasture, Agriculture, Development, Water, cover_2022)
gc()


# 2021 ----
Agriculture <- cover_2021
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- FALSE
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- TRUE
Agriculture <- raster(Agriculture)

Development <- cover_2021
Development[Development %notin% c(24,25,30)] <- FALSE
Development[Development %in% c(24,25,30)] <- TRUE
Development <- raster(Development)

Forest <- cover_2021
Forest[Forest %notin% c(1,3,4,5,6,49,29,9)] <- FALSE
Forest[Forest %in% c(1,3,4,5,6,49,29,9)] <- TRUE
Forest <- raster(Forest)

Pasture <- cover_2021
Pasture[Pasture %notin% c(12,15)] <- FALSE
Pasture[Pasture %in% c(12,15)] <- TRUE
Pasture <- raster(Pasture)

Water <- cover_2021
Water[Water %notin% c(11,26,33)] <- FALSE
Water[Water %in% c(11,26,33)] <- TRUE
Water <- raster(Water)

covers_2021 <- list(Agriculture = Agriculture,
                    Development = Development,
                    Forest = Forest,
                    Pasture = Pasture,
                    Water = Water)

# free environment space
rm(Native_Forest, Forestry, Pasture, Agriculture, Development, Water, cover_2021)
gc()


# 2020 ----
Agriculture <- cover_2020
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- FALSE
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- TRUE
Agriculture <- raster(Agriculture)

Development <- cover_2020
Development[Development %notin% c(24,25,30)] <- FALSE
Development[Development %in% c(24,25,30)] <- TRUE
Development <- raster(Development)

Forest <- cover_2020
Forest[Forest %notin% c(1,3,4,5,6,49,29,9)] <- FALSE
Forest[Forest %in% c(1,3,4,5,6,49,29,9)] <- TRUE
Forest <- raster(Forest)

Pasture <- cover_2020
Pasture[Pasture %notin% c(12,15)] <- FALSE
Pasture[Pasture %in% c(12,15)] <- TRUE
Pasture <- raster(Pasture)

Water <- cover_2020
Water[Water %notin% c(11,26,33)] <- FALSE
Water[Water %in% c(11,26,33)] <- TRUE
Water <- raster(Water)

covers_2020 <- list(Agriculture = Agriculture,
                    Development = Development,
                    Forest = Forest,
                    Pasture = Pasture,
                    Water = Water)

# free environment space
rm(Native_Forest, Forestry, Pasture, Agriculture, Development, Water, cover_2020)
gc()


# 2019 ----
Agriculture <- cover_2019
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- FALSE
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- TRUE
Agriculture <- raster(Agriculture)

Development <- cover_2019
Development[Development %notin% c(24,25,30)] <- FALSE
Development[Development %in% c(24,25,30)] <- TRUE
Development <- raster(Development)

Forest <- cover_2019
Forest[Forest %notin% c(1,3,4,5,6,49,29,9)] <- FALSE
Forest[Forest %in% c(1,3,4,5,6,49,29,9)] <- TRUE
Forest <- raster(Forest)

Pasture <- cover_2019
Pasture[Pasture %notin% c(12,15)] <- FALSE
Pasture[Pasture %in% c(12,15)] <- TRUE
Pasture <- raster(Pasture)

Water <- cover_2019
Water[Water %notin% c(11,26,33)] <- FALSE
Water[Water %in% c(11,26,33)] <- TRUE
Water <- raster(Water)

covers_2019 <- list(Agriculture = Agriculture,
                    Development = Development,
                    Forest = Forest,
                    Pasture = Pasture,
                    Water = Water)

# free environment space
rm(Native_Forest, Forestry, Pasture, Agriculture, Development, Water, cover_2019)
gc()


# 2018 ----
Agriculture <- cover_2018
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- FALSE
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- TRUE
Agriculture <- raster(Agriculture)

Development <- cover_2018
Development[Development %notin% c(24,25,30)] <- FALSE
Development[Development %in% c(24,25,30)] <- TRUE
Development <- raster(Development)

Forest <- cover_2018
Forest[Forest %notin% c(1,3,4,5,6,49,29,9)] <- FALSE
Forest[Forest %in% c(1,3,4,5,6,49,29,9)] <- TRUE
Forest <- raster(Forest)

Pasture <- cover_2018
Pasture[Pasture %notin% c(12,15)] <- FALSE
Pasture[Pasture %in% c(12,15)] <- TRUE
Pasture <- raster(Pasture)

Water <- cover_2018
Water[Water %notin% c(11,26,33)] <- FALSE
Water[Water %in% c(11,26,33)] <- TRUE
Water <- raster(Water)

covers_2018 <- list(Agriculture = Agriculture,
                    Development = Development,
                    Forest = Forest,
                    Pasture = Pasture,
                    Water = Water)

# free environment space
rm(Native_Forest, Forestry, Pasture, Agriculture, Development, Water, cover_2018)
gc()


# 2017 ----
Agriculture <- cover_2017
Agriculture[Agriculture %notin% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- FALSE
Agriculture[Agriculture %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48)] <- TRUE
Agriculture <- raster(Agriculture)

Development <- cover_2017
Development[Development %notin% c(24,25,30)] <- FALSE
Development[Development %in% c(24,25,30)] <- TRUE
Development <- raster(Development)

Forest <- cover_2017
Forest[Forest %notin% c(1,3,4,5,6,49,29,9)] <- FALSE
Forest[Forest %in% c(1,3,4,5,6,49,29,9)] <- TRUE
Forest <- raster(Forest)

Pasture <- cover_2017
Pasture[Pasture %notin% c(12,15)] <- FALSE
Pasture[Pasture %in% c(12,15)] <- TRUE
Pasture <- raster(Pasture)

Water <- cover_2017
Water[Water %notin% c(11,26,33)] <- FALSE
Water[Water %in% c(11,26,33)] <- TRUE
Water <- raster(Water)

covers_2017 <- list(Agriculture = Agriculture,
                    Development = Development,
                    Forest = Forest,
                    Pasture = Pasture,
                    Water = Water)

# free environment space
rm(Native_Forest, Forestry, Pasture, Agriculture, Development, Water, cover_2017)
gc()

