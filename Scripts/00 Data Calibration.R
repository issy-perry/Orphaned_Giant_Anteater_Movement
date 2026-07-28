# 0 Data Calibration Script
# This script includes the following:
#     1. Extracting telemetry points from non-moving points
#     2. Estimate UERE from non-moving points
#     3. Apply new UERE estimates to all other moving points
# now we can use the telemetry points from known mortalities as calibration points.
# we need to do this because the collar was likely calibrated on an individual within a pen (little bits of movement can still severly impact the calibration of these devices)
# this will result with more accurately calibrated data, as we know that the collar was not moving
# wild-raised individuals' data were pre-calibrated as it was used in another project

# load packages
library(ctmm) # converting data into a telemetry object
library(tidyverse) # formatting of tbl for Stefano's functions
library(sf)
library(lubridate) # included in the structure of Stefano's functions
library(mapview) # included in the structure of Stefano's functions


# reimport data for individuals who died
Bahia_gps <- read.csv ("./DATA/Raw_Data/New_Orphans/Bahia/729088A_1_Datalog.csv", skip = 22)
Capitu_gps <- read.csv ("./DATA/Raw_Data/Old_Orphans/709201A_Datalog.csv", skip = 22)
Cláudio_gps <- read.csv ("./DATA/Raw_Data/New_Orphans/Cláudio/716785A_2_Datalog.csv", skip = 22)
Mulan_gps <- read.csv ("./DATA/Raw_Data/New_Orphans/Mulan/716786A_1_Datalog.csv", skip = 22)
Nayeli_gps <- read.csv ("./DATA/Raw_Data/New_Orphans/Nayeli/729091A_1_Datalog.csv", skip = 22)
RandN_gps <-read.csv ("./DATA/Raw_Data/New_Orphans/Rita_e_Nancy/729101A_1_Datalog.csv", skip = 22)
# skipping Arya, as she drowned in a pool and may not have been completely stagnate

# remove NA values using the GPS Latitude column (no successful fix should have an NA for location)
Bahia_gps <- Bahia_gps [!is.na(Bahia_gps$GPS.Latitude),]
Capitu_gps <- Capitu_gps [!is.na(Capitu_gps$GPS.Latitude),]
Cláudio_gps <- Cláudio_gps [!is.na(Cláudio_gps$GPS.Latitude),]
Mulan_gps <- Mulan_gps [!is.na(Mulan_gps$GPS.Latitude),]
Nayeli_gps <- Nayeli_gps [!is.na(Nayeli_gps$GPS.Latitude),]
RandN_gps <- RandN_gps [!is.na(RandN_gps$GPS.Latitude),]

# reformat time stamp
Bahia_gps$GPS.Fix.Time <- as.POSIXct(Bahia_gps$GPS.Fix.Time, format = "%Y.%m.%d %H:%M:%S")
Capitu_gps$GPS.Fix.Time <- as.POSIXct(Capitu_gps$GPS.Fix.Time, format = "%Y.%m.%d %H:%M:%S")
Cláudio_gps$GPS.Fix.Time <- as.POSIXct(Cláudio_gps$GPS.Fix.Time, format = "%Y.%m.%d %H:%M:%S")
Mulan_gps$GPS.Fix.Time <- as.POSIXct(Mulan_gps$GPS.Fix.Time, format = "%Y.%m.%d %H:%M:%S")
Nayeli_gps$GPS.Fix.Time <- as.POSIXct(Nayeli_gps$GPS.Fix.Time, format = "%Y.%m.%d %H:%M:%S")
RandN_gps$GPS.Fix.Time <- as.POSIXct(RandN_gps$GPS.Fix.Time, format = "%Y.%m.%d %H:%M:%S")

# it is impossible these anteaters to be above 0 latitudinally, so we can safely exclude these points
Bahia_gps <- Bahia_gps[Bahia_gps$GPS.Latitude < 0,]
Capitu_gps <- Capitu_gps[Capitu_gps$GPS.Latitude < 0,]
Cláudio_gps <- Cláudio_gps[Cláudio_gps$GPS.Latitude < 0,]
Mulan_gps <- Mulan_gps[Mulan_gps$GPS.Latitude < 0,]
Nayeli_gps <- Nayeli_gps[Nayeli_gps$GPS.Latitude < 0,]
RandN_gps <- RandN_gps[RandN_gps$GPS.Latitude < 0,]

# only Rita's data will be used, so her points need to be subset to not include Nancy's
Rita_gps <- RandN_gps [RandN_gps$GPS.Fix.Time >= ("2022-12-14 13:20:26") & RandN_gps$GPS.Fix.Time <= ("2023-02-12 16:01:31"),]
rm(RandN_gps)

# remove predeployment points
Bahia_gps <- Bahia_gps[Bahia_gps$Predeployment.Data == "No",]
Capitu_gps <- Capitu_gps[Capitu_gps$Predeployment.Data == "No",]
Cláudio_gps <- Cláudio_gps[Cláudio_gps$Predeployment.Data == "No",]
Mulan_gps <- Mulan_gps[Mulan_gps$Predeployment.Data == "No",]
Nayeli_gps <- Nayeli_gps[Nayeli_gps$Predeployment.Data == "No",]
Rita_gps <- Rita_gps[Rita_gps$Predeployment.Data == "No",]


# only include mortality data
Bahia_gps <- Bahia_gps[Bahia_gps$Mortality == "Yes*",]
Capitu_gps <- Capitu_gps[Capitu_gps$Mortality == "Yes*",]
Cláudio_gps <- Cláudio_gps[Cláudio_gps$Mortality == "Yes*",]
Mulan_gps <- Mulan_gps[Mulan_gps$Mortality == "Yes*",]
Nayeli_gps <- Nayeli_gps[Nayeli_gps$Mortality == "Yes*",]
Rita_gps <- Rita_gps[Rita_gps$Mortality == "Yes*",]

# remove if there was a clear movement/translocation
Bahia_gps <- Bahia_gps[Bahia_gps$GPS.Fix.Time <= ("2024-02-24 18:40:16"),] 

# drop Cláudio, as he does not have enough fixes (only 3)
rm(Cláudio_gps)

# make an outlier column
Bahia_gps$outlier <- replicate (59, 0)
Capitu_gps$outlier <- replicate(48, 0)
Mulan_gps$outlier <- replicate(103, 0)
Nayeli_gps$outlier <- replicate (23, 0)
Rita_gps$outlier <- replicate (19, 0)

# make another tibble to see if the animals are being transported or moving, as any movements would disrupt calibration
d <- tibble(animal = c('Bahia', 'Capitu', 'Mulan', 'Nayeli', 'Rita'),
            tel = purrr::map(animal, function(.name) get(paste0(.name, '_gps'))))

# check each individual to make sure their bodies were not moved (we only want to see one cluster of points with some outlier spiking)
# Bahia
out <- check_animal("Bahia") # just outliers

# Capitu
out <- check_animal("Capitu") # just outliers

# Mulan
out <- check_animal("Mulan") # just outliers

# Nayeli
out <- check_animal("Nayeli") # just outliers

# Rita
out <- check_animal("Rita") # just outliers

# save mortality output
saveRDS(object = d,
        file = paste0("./DATA/mortality_data",
                      format(Sys.time(), '%Y-%m-%d-%H-%M'),
                      '.rds'))



# import mortality data and convert to telemetry 
mortality_tbl <- readRDS("./DATA/mortality_data2025-04-28-20-39.rds")

# check class
class(mortality_tbl) #dataset was imported as a tbl_df

# add IDs back in
for(i in 1:length(mortality_tbl$tel)){ 
  mortality_tbl$tel[[i]]$ID <- mortality_tbl$animal[[i]] 
}

# convert the tbl into a dataframe
mortality <- do.call(rbind, mortality_tbl$tel) #converts tiblle to a dataframe

# check class to see if change wen through
class(mortality) #works

# convert the dataframe to a ctmm telemetry object (will output a nested list)
mortality <- as.telemetry(mortality, mark.rm = TRUE) # convert to telemetry object and list and mark.rm = TRUE will drop outliers

# save file
save(mortality, file = "./DATA/Orphaned_tel_data/Mortality_tel.rda") 

# estimate UERE from mortality points ----
# load mortality and orphan data
load("./DATA/Orphaned_tel_data/Data_telemetry_NO_UERE.rda") #orphaned telemetry data
load("./DATA/Orphaned_tel_data/Mortality_tel.rda") #mortality telemetry data

# estimate the RMS UERE for the dead individuals
UERE_2 <- uere.fit(mortality)

# plot example to see difference
plot(DATA_orphan[["Dom"]], main = "Before")

# assign the estimated UERE to the whole dataset
uere(DATA_orphan) <- UERE_2

# plot difference
plot(DATA_orphan[["Dom"]], main = "After")

# save telemetry data
save(DATA_orphan, file = "./DATA/Orphaned_tel_data/Data_telemetry.rda") 