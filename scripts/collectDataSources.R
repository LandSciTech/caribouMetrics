# Best available data required to assess the observed state of the Churchill
# range for caribou in 2010 and 2020. There are discrepancies among data sources
# describing the state of the range, so ECCC will provide up to 4 different
# descriptions of the current and past state of the range.

# Data required:
## Static:
### PLC
### Eskers

## Dynamic:
### FRI/Age 2010 and 2020
### Natural disturbance ie Fire 2010 and 2020
### Harvest disturbance 2010 and 2020
### Roads, Rail and Utilities 2010 and 2020

# Data Sources
friD <- st_read("inputNV/ChurchillFRI/FMU_120trlake_175caribou_702lacseul/FMU_120trlake_175caribou_702lacseul.shp")
# Seems to be from 2012 based on Age + YrOrg Seems to be most recent available
# on GeoHub but check with Josie

# Roads
# From Rob Rempel
roadRempel <- st_read("inputNV/HornsethRempelInfo/Dist_UtilLinesEtc_Churchill/Roads_Churchill.shp")
# most recent year construtructed 2009

# From Josie
road2 <- st_read("inputNV/Roads/roads.gdb", layer = "MNRFwithORN")
# most recent year is 2017

# From https://geohub.lio.gov.on.ca/datasets/mnrf-road-segments and
# https://geohub.lio.gov.on.ca/datasets/mnrf::ontario-road-network-orn-segment-with-address

# all of ON up to 2020
road_ORN <- st_read("inputNV/Roads/ORNSEGAD_20201028/Non_Sensitive.gdb")

caribouRanges <- st_read("./inputNV/caribouRanges/Caribou_Range_Boundary.shp", 
                        quiet = TRUE) %>% 
  st_transform(st_crs(road_ORN))

road_MNRF <- st_read("inputNV/Roads/MNRRDSEG_20201028/Non_Sensitive.gdb", 
                     layer = "MNRF_ROAD_SEGMENT")
road_MNRF2 <- st_filter(road_MNRF, caribouRanges)

road_ORN2 <- st_filter(road_ORN, caribouRanges)

road_MNRF_ORN <- rbind(road_MNRF2, road_ORN2)

# Utilities use Business effective date for age "Date that the record becomes
# effective in relation to the business i.e. the date MNR became aware of its
# existence."
utilRempel <- st_read("inputNV/HornsethRempelInfo/Dist_UtilLinesEtc_Churchill/UtilLines_Churchill.shp")
# last Business_effective_date is 2007

# from GeoHub
utilON <- st_read("inputNV/Utilities/Utility_Line.shp")
#last business effective date is 2019 colname is BUSINESS_1
# Note some utilities are overlaping utilities of different types

utilONCH <- utilON %>% st_crop(projectPolyD) 

# the ones in Rempel and not ONCH are all on the edge
plot(st_bbox(projectPolyD) %>% st_as_sfc() %>% st_transform(st_crs(utilRempel)))
anti_join(utilRempel,
          utilONCH %>% st_drop_geometry(), by = "OGF_ID") %>%
  select(OGF_ID) %>% plot(add = T)

railRempel <- st_read("inputNV/HornsethRempelInfo/Dist_UtilLinesEtc_Churchill/Railyway_Churchill.shp")

railON <- st_read("inputNV/Rail/ORWN_TRACK.shp")

railONCH <- railON %>% st_transform(st_crs(utilON)) %>% st_crop(projectPolyD) 

# no dates with railways but have new ones been built?

#Disturbance #================================================================
pthBase <- "inputNV/Disturbance/"

# BEAD 
## Only available 2010 and 2015 as snap shot of currently visible disturbance
## plus fires for previous 40 years
devtools::load_all(".")
# NRCAN fire data base
firesNFDB <- st_read(paste0(pthBase,
                          "NFDB_poly/NFDB_poly_20201005.shp"))

caribouRanges <- st_read("./inputNV/caribouRanges/Caribou_Range_Boundary.shp", 
                         quiet = TRUE)

firesNFDBCaribou <- firesNFDB %>% 
  st_make_valid() %>%
  st_filter(caribouRanges %>% st_transform(crs = st_crs(firesNFDB)))


st_write(firesNFDBCaribou,
         paste0(pthBase, "NFDB_poly/NFDB_ON_Caribou.shp"),
         append = FALSE)

firesNFDBCaribou <- st_read( paste0(pthBase, "NFDB_poly/NFDB_ON_Caribou.shp"))

# CanLaD Canada Landsat Disturbance
## NRCAN data that covers 1984 to 2015 fire and harvest
## https://open.canada.ca/data/en/dataset/add1346b-f632-4eb9-a83d-a662b38655ad

CanLadYr <- raster("inputNV/Disturbance/CanLaD/CanLaD_20151984_latest_YRT2.tif")
CanLadYrCaribou <- raster::crop(CanLadYr, 
                                caribouRanges %>% st_transform(st_crs(CanLadYr)),
                                filename = "inputNV/Disturbance/CanLaD/Caribou_CanLaDYr.tif")

# High resolution forest change for Canada
## NRCAN data that covers 1985 to 2015 fire and harvest
## https://open.canada.ca/data/en/dataset/eaf92ab8-c744-4f8f-a280-6987d4f6005e
## this seems to be the data used to create the OFRID Landsat data for 1984-2015

hRFCCanYr <- raster("inputNV/Disturbance/C2C_change_year_2012_2015/C2C_change_year_2012_2015.tif")

hRFCONCaribouYr <- raster::crop(hRFCCanYr, 
                              caribouRanges %>% st_transform(st_crs(hRFCCanYr)),
                              filename = "inputNV/Disturbance/C2C_change_year_2012_2015/C2C_change_year_2012_2015_ONCaribou.tif")

hRFCCanTy <- raster("inputNV/Disturbance/C2C_change_type_2012_2015/C2C_change_type_2012_2015.tif")

hRFCONCaribouTy <- raster::crop(hRFCCanTy, 
                              caribouRanges %>% st_transform(st_crs(hRFCCanTy)),
                              filename = "inputNV/Disturbance/C2C_change_type_2012_2015/C2C_change_type_2012_2015_ONCaribou.tif",
                              overwrite = TRUE)
## type 0 = no change, 1 = wildfire, 2 = harvest, 5 = lower confidence wildfire,
## 6 = lower confidence harvest

plot(hRFCONCaribouTy == 1|hRFCONCaribouTy == 5)

plot(firesNFDBCaribou %>% filter(YEAR >= 2012, YEAR <= 2015) %>% st_geometry(), add = T)

# OFRI Disturbance database Use original data sources instead
# MNRF/AFFES https://geohub.lio.gov.on.ca/datasets/fire-disturbance-area
AFFES <- st_read("inputNV/Disturbance/FIREDSTB/LIO-2020-10-28/FIRE_DISTURBANCE_AREA.shp")


# Far North land cover
## Newer land cover data doesn't seem to fully overlay churchill but does for
## the rest. Similar classes as PLC2000 but supposed to be more accurate.
## However, diferent classification methods so wouldn't have same biases.
fNLC <- raster("inputNV/Provincial-Landcover-2000/FarNorthLandCover/Version 1.4/TIF Format/Class/FarNorth_LandCover_Class_UTM16.tif")

fNLCCaribou <- raster::crop(fNLC, 
                                caribouRanges %>% st_transform(st_crs(fNLC)),
                                filename = "inputNV/Provincial-Landcover-2000/FarNorthLandCover/UTM15_Caribou.tif",
                                overwrite = TRUE)
