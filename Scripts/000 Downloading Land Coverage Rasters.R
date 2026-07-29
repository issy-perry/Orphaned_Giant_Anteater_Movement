# 000 Downloading Land Coverage Rasters
# This script includes the following:
#     1. Download rasters for entire country of Brazil
#     2. Extract extents of telemetry data for cropping the rasters
#     3. Crop rasters


# load packages
library(terra) # main package we need to use for dowloading, cropping, and saving our rasters
library(sf) # getting the extents of datasets

# load in land coverage data
land_cover_2024 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2024.tif")
land_cover_2023 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2023.tif")
land_cover_2022 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2022.tif")
land_cover_2021 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2021.tif")
land_cover_2020 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2020.tif")
land_cover_2019 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2019.tif")
land_cover_2018 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2018.tif")
land_cover_2017 <- terra::rast("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_10/lulc/coverage/brazil_coverage_2017.tif")

# load tel data
load("./DATA/Orphaned_tel_data/Data_telemetry_RR.rda") # orphaned 
load("./DATA/Wild_raised_tel_data/Data_telemetry_RR.rda") # wild-raised

# get extents from tel data to use for cropping the rasters
# wild-raised
# convert back to a dataframe
DATA_wild <- do.call(rbind.data.frame, DATA_wild)
# convert to a sf object to get extents of tel data
DATA_wild_sf <- st_as_sf(DATA_wild, coords = c("longitude", "latitude"), 
                         crs = crs(land_cover_2017))
# get extent barrier of tel data
ew <- ext(DATA_wild_sf)
# convert to a polygon 
POLY <- as.polygons(DATA_wild_ext, crs = "EPSG:4326") 
# create a 4000 m buffer around the polygon (buffer ensures all points and UD boundaries are within the raster's boundary)
POLY_buff <- buffer(POLY, width = 4000)
# create new extent that includes buffer
ew <- ext(POLY_buff)
# print to check
print(ew)


# orphaned
# convert back to a dataframe
DATA_orphan <- do.call(rbind.data.frame, DATA_orphan)
# convert to an sf object to get extents of tel data
DATA_orphan_sf <- st_as_sf(DATA_orphan, coords = c("longitude", "latitude"), 
                         crs = crs(land_cover_2019))
# get extent barrier of tel data
eo <- ext(DATA_orphan_sf)
# convert to a polygon 
POLY <- as.polygons(DATA_orphan_ext, crs = "EPSG:4326") 
# create a 4000 m buffer around the polygon (buffer ensures all points and UD boundaries are within the raster's boundary)
POLY_buff <- buffer(POLY, width = 4000)
# create new extent that includes buffer
eo <- ext(POLY_buff)
# print to check
print(eo)

# crop rasters using previously made extents
cover_2024 <- crop(land_cover_2024, eo)
cover_2023 <- crop(land_cover_2023, eo)
cover_2022 <- crop(land_cover_2022, eo)
cover_2021 <- crop(land_cover_2021, eo)
cover_2020 <- crop(land_cover_2020, eo)
cover_2019 <- crop(land_cover_2019, eo)
cover_2018 <- crop(land_cover_2018, ew)
cover_2017 <- crop(land_cover_2017, ew)

# free environment space
rm(land_cover_2024, land_cover_2023, land_cover_2022, land_cover_2021, land_cover_2020, land_cover_2019, land_cover_2018, land_cover_2017, ew, eo)
gc()



