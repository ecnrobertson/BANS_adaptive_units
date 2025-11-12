library(terra)

BANS.trend <- rast("data/spatial_files/Trend maps raster tif files/tr06160.tif")
plot(BANS.trend)

res(BANS.trend)


library(openxlsx)

master <- read.xlsx("~/Desktop/GCRF/Zavaleta_Brown_SNRF_Project_BandingData_2018-2022.xlsx")

library(dplyr)
actual_samples <- read.csv("~/OneDrive - Colostate/GCRF_Pub/data/raw_GCRF_data_all_morpho.csv")
TB <- master %>% filter(band_number %in% actual_samples$band_number) %>% filter(bander=="TB")
nrow(TB)
EZ <- master %>% filter(band_number %in% actual_samples$band_number) %>% filter(bander=="EZ")
nrow(EZ)
