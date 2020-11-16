# create data sets for Churchill and ROF

devtools::load_all(".")

outCHill <- "inputNV/ChurchillData/"
outROF <- "inputNV/ROFData/"
outBSand <- "inputNV/BrightsandData/"

caribouRanges <- st_read("./inputNV/caribouRanges/Caribou_Range_Boundary.shp", 
                         quiet = TRUE)

# use crs from eskers
eskerD <- st_read("inputNV/HornsethRempelInfo/Eskers_Ontario/Eskers_Ontario.shp")

crsUse <- st_crs(eskerD) 
crsUseR <- crs(eskerD)
# buffer project areas to 6000m because that is larger than the largest window
# radius times 3 to account for offsets

rof <- caribouRanges %>% 
  filter(RANGE_NAME %in% c("Missisa", "James Bay", 
                           "Pagwachuan", "Nipigon", "Ozhiski")) %>% 
  st_transform(crsUse) %>% 
  summarise(OGF_ID = first(OGF_ID),
            geometry = st_union(geometry) %>% st_buffer(6000 *3))


cHill <- caribouRanges %>% filter(RANGE_NAME == "Churchill") %>% 
  st_transform(crsUse) %>% 
  transmute(OGF_ID = first(OGF_ID),
            geometry = st_buffer( geometry, 6000 *3))

bSand <- caribouRanges %>% filter(RANGE_NAME == "Brightsand") %>% 
  st_transform(crsUse) %>% 
  transmute(OGF_ID = first(OGF_ID),
            geometry = st_buffer( geometry, 6000 *3))

# plc #=========================================================================
# No change use the same for all

# 250 or 15? Only have 15 for whole province
# plc15 <- raster("inputNV/Provincial-Landcover-2000/Provincial Landcover 2000/z15-27class.tif")
# plc16 <- raster("inputNV/Provincial-Landcover-2000/Provincial Landcover 2000/z16-27class.tif") 
# plc17 <- raster("inputNV/Provincial-Landcover-2000/Provincial Landcover 2000/z17-27class.tif") 
# 
# tmplt15 <- raster::projectExtent(plc15, crs = crs(plc16)) %>% raster::extent()
# tmplt16 <- raster::extent(plc16)
# tmplt17 <- raster::projectExtent(plc17, crs = crs(plc16)) %>% raster::extent()
# 
# 
# tmplt <- raster(raster::extent(tmplt15@xmin, tmplt17@xmax, 
#                 min(tmplt15@ymin, tmplt16@ymin, tmplt17@ymin),
#                 max(tmplt15@ymax, tmplt16@ymax, tmplt17@ymax))) %>% 
#   raster::`crs<-`(value = crs(plc16)) %>% 
#   raster::`res<-`(value = res(plc16)) %>% 
#   raster::projectExtent(crsUseR)

# merged in QGIS 

plcON <- raster("inputNV/Provincial-Landcover-2000/Provincial Landcover 2000/z15-17merged.tif")

plcCH <- raster::crop(plcON, cHill %>% st_transform(st_crs(plcON)), 
                      filename = paste0(outCHill, "plc.tif"),
                      datatype = "INT1U")

plcROF <- raster::crop(plcON, rof %>% st_transform(st_crs(plcON)), 
                      filename = paste0(outROF, "plc.tif"),
                      datatype = "INT1U")

plcBSand <- raster::crop(plcON, bSand %>% st_transform(st_crs(plcON)), 
                         filename = paste0(outBSand, "plc.tif"),
                         datatype = "INT1U")

# Far North land cover is available but not currently using

## Newer land cover data doesn't seem to fully overlay churchill but does for
## the rest. Similar classes as PLC2000 but supposed to be more accurate.
## However, diferent classification methods so wouldn't have same biases.

# fNLC <- raster("inputNV/Provincial-Landcover-2000/FarNorthLandCover/Version 1.4/TIF Format/Class/FarNorth_LandCover_Class_UTM16.tif")
# 
# fNLCCaribou <- raster::crop(fNLC, 
#                             caribouRanges %>% st_transform(st_crs(fNLC)),
#                             filename = "inputNV/Provincial-Landcover-2000/FarNorthLandCover/UTM15_Caribou.tif",
#                             overwrite = TRUE)

# eskers #======================================================================
eskerD <- st_read("inputNV/HornsethRempelInfo/Eskers_Ontario/Eskers_Ontario.shp")

eskerCH <- st_filter(eskerD, cHill)

eskerROF <- st_filter(eskerD, rof)

eskerBSand <- st_filter(eskerD, bSand)

st_write(eskerCH, paste0(outCHill, "esker.shp"))
st_write(eskerROF, paste0(outROF, "esker.shp"))
st_write(eskerBSand, paste0(outBSand, "esker.shp"))

rm(eskerD, eskerCH, eskerROF, eskerBSand)

# road #========================================================================

# Need 2010 and 2020 versions use year constructed in MNRF data
road_ORN <- st_read("inputNV/Roads/ORNSEGAD_20201028/Non_Sensitive.gdb")

road_MNRF <- st_read("inputNV/Roads/MNRRDSEG_20201028/Non_Sensitive.gdb", 
                     layer = "MNRF_ROAD_SEGMENT")

# Churchill
road_ORNCH <- st_filter(road_ORN, cHill %>% st_transform(st_crs(road_ORN))) %>% 
  st_transform(crsUse)

road_MNRFCH <- st_filter(road_MNRF, cHill %>% st_transform(st_crs(road_MNRF))) %>% 
  st_transform(crsUse)

st_write(road_MNRFCH %>% st_cast("LINESTRING"),
         paste0(outCHill, "temproad_MNRFCH.gpkg"), append = FALSE)
st_write(road_ORNCH %>% st_cast("LINESTRING"), 
         paste0(outCHill, "temproad_ORNCH.gpkg"), append = FALSE)

# Used QGIS to buffer MNRF by 100 m and then difference that with ORN. So ORN -
# MNRF, then merge them together and saved as road_ORNMNRF.gpkg

road_ORNMNRFCH <- st_read(paste0(outCHill, "road_ORNMNRFCH.gpkg")) %>% 
  st_set_crs(crsUse) %>% select(YEAR_CONSTRUCTED)

road_ORNMNRFCH2010 <- road_ORNMNRFCH %>% 
  filter(YEAR_CONSTRUCTED <= 2010 | is.na(YEAR_CONSTRUCTED))

st_write(road_ORNMNRFCH, paste0(outCHill, "road_ORNMNRFCH2020.shp"))
st_write(road_ORNMNRFCH2010, paste0(outCHill, "road_ORNMNRFCH2010.shp"))

road_ORNMNRFCH2020 <- st_read(paste0(outCHill, "road_ORNMNRFCH2020.shp"))
road_ORNMNRFCH2010 <- st_read(paste0(outCHill, "road_ORNMNRFCH2010.shp"))

plot(road_ORNMNRFCH2020 %>% st_geometry(), col = "red")
plot(road_ORNMNRFCH2010 %>% st_geometry(), add = TRUE)
plot(road_ORNMNRFCH2020 %>% filter(is.na(YEAR_CO)) %>% st_geometry(), 
     col = "green", add = TRUE)
plot(road_ORNMNRFCH2020 %>% filter(YEAR_CO == 9999) %>% st_geometry(), 
     col = "blue", add = TRUE)

rm(road_ORNMNRFCH, road_MNRFCH, road_ORNCH, road_ORNMNRFCH2010, 
   road_ORNMNRFCH2020)

# ROF
road_ORNROF <- st_filter(road_ORN, rof %>% st_transform(st_crs(road_ORN))) %>% 
  st_transform(crsUse)

road_MNRFROF <- st_filter(road_MNRF, rof %>% st_transform(st_crs(road_MNRF))) %>% 
  st_transform(crsUse)

st_write(road_MNRFROF %>% st_cast("LINESTRING"),
         paste0(outROF, "temproad_MNRFROF.gpkg"), append = FALSE)
st_write(road_ORNROF %>% st_cast("LINESTRING"), 
         paste0(outROF, "temproad_ORNROF.gpkg"), append = FALSE)

# Used QGIS to buffer MNRF by 100 m and then difference that with ORN. So ORN -
# MNRF, then merge them together and saved as road_ORNMNRFROF.gpkg

road_ORNMNRFROF2020 <- st_read(paste0(outROF, "road_ORNMNRFROF.gpkg")) %>% 
  st_set_crs(crsUse) %>% select(YEAR_CONSTRUCTED)

road_ORNMNRFROF2010 <- road_ORNMNRFROF2020 %>% 
  filter(YEAR_CONSTRUCTED <= 2010 | is.na(YEAR_CONSTRUCTED))

st_write(road_ORNMNRFROF2020, paste0(outROF, "road_ORNMNRFROF2020.shp"))
st_write(road_ORNMNRFROF2010, paste0(outROF, "road_ORNMNRFROF2010.shp"))

plot(road_ORNMNRFROF2020 %>% st_geometry(), col = "red")
plot(road_ORNMNRFROF2010 %>% st_geometry(), add = TRUE)
plot(road_ORNMNRFROF2020 %>% filter(is.na(YEAR_CONSTRUCTED)) %>% st_geometry(), 
     col = "green", add = TRUE)
plot(road_ORNMNRFROF2020 %>% filter(YEAR_CONSTRUCTED == 9999) %>% st_geometry(), 
     col = "blue", add = TRUE)

# Brightsand
road_ORNBSand <- st_filter(road_ORN, bSand %>% st_transform(st_crs(road_ORN))) %>% 
  st_transform(crsUse)

road_MNRFBSand <- st_filter(road_MNRF, bSand %>% st_transform(st_crs(road_MNRF))) %>% 
  st_transform(crsUse)

st_write(road_MNRFBSand %>% st_cast("LINESTRING"),
         paste0(outBSand, "temproad_MNRFBSand.gpkg"), append = FALSE)
st_write(road_ORNBSand %>% st_cast("LINESTRING"), 
         paste0(outBSand, "temproad_ORNBSand.gpkg"), append = FALSE)

# Used QGIS to buffer MNRF by 100 m and then difference that with ORN. So ORN -
# MNRF, then merge them together and saved as road_ORNMNRFROF.gpkg

road_ORNMNRFBSand2020 <- st_read(paste0(outBSand, "road_ORNMNRFBSand.gpkg")) %>% 
  st_set_crs(crsUse) %>% select(YEAR_CONSTRUCTED)

road_ORNMNRFBSand2010 <- road_ORNMNRFBSand2020 %>% 
  filter(YEAR_CONSTRUCTED <= 2010 | is.na(YEAR_CONSTRUCTED))

st_write(road_ORNMNRFBSand2020, paste0(outBSand, "road_ORNMNRFBSand2020.shp"))
st_write(road_ORNMNRFBSand2010, paste0(outBSand, "road_ORNMNRFBSand2010.shp"))

plot(road_ORNMNRFBSand2020 %>% st_geometry(), col = "red")
plot(road_ORNMNRFBSand2010 %>% st_geometry(), add = TRUE)
plot(road_ORNMNRFBSand2020 %>% filter(is.na(YEAR_CONSTRUCTED)) %>% st_geometry(), 
     col = "green", add = TRUE)
plot(road_ORNMNRFBSand2020 %>% filter(YEAR_CONSTRUCTED == 9999) %>% st_geometry(), 
     col = "blue", add = TRUE)

rm(road_MNRF, road_MNRFBSand, road_MNRFROF, road_ORN, road_ORNBSand, 
   road_ORNMNRFBSand2010, road_ORNMNRFBSand2020, road_ORNMNRFROF2010,
   road_ORNMNRFROF2020, road_ORNROF)

# rail #========================================================================
railON <- st_read("inputNV/Rail/ORWN_TRACK.shp") %>% 
  st_transform(crsUse)

# rail does not have dates associated so assuming no new ones
# seems correct for Churchill when compared to what Rempel used

railCH <- st_filter(railON, cHill)

railROF <- st_filter(railON, rof)

railBSand <- st_filter(railON, bSand)

st_write(railCH, paste0(outCHill, "rail.shp"))
st_write(railROF, paste0(outROF, "rail.shp"))
st_write(railBSand, paste0(outBSand, "rail.shp"))

rm(railON, railCH, railROF, railBSand)

# Utilities #===================================================================
utilON <- st_read("inputNV/Utilities/Utility_Line.shp") %>% 
  st_transform(crsUse)

# Utilities use Business effective date for age "Date that the record becomes
# effective in relation to the business i.e. the date MNR became aware of its
# existence."

# Note some utilities are partially overlaping utilities of different types 
# not clear what Rempel did about this. Also often follow roads

utilCH2020 <- st_filter(utilON, cHill)

utilCH2010 <- utilCH2020 %>% filter(BUSINESS_1 < as.Date("2010-01-01"))

utilROF2020 <- st_filter(utilON, rof)

utilROF2010 <- utilROF2020 %>% filter(BUSINESS_1 < as.Date("2010-01-01"))

utilBSand2020 <- st_filter(utilON, bSand)

utilBSand2010 <- utilBSand2020 %>% filter(BUSINESS_1 < as.Date("2010-01-01"))

plot(utilROF2020 %>% st_geometry())
plot(utilROF2010 %>% st_geometry(), add = TRUE, col = "red")

st_write(utilCH2020, paste0(outCHill, "util2020.shp"))
st_write(utilROF2020, paste0(outROF, "util2020.shp"))
st_write(utilBSand2020, paste0(outBSand, "util2020.shp"))

st_write(utilCH2010, paste0(outCHill, "util2010.shp"))
st_write(utilROF2010, paste0(outROF, "util2010.shp"))
st_write(utilBSand2010, paste0(outBSand, "util2010.shp"))

rm(utilON, utilCH2020, utilROF2020, utilBSand2020, 
   utilCH2010, utilROF2010, utilBSand2010)

# fri #=======================================================================
# Waiting on updated FRI data from Benoit

# Have Churchill Not sure when it is current to age + year origin = 2012 
# based on https://www.sdc.gov.on.ca/sites/MNRF-PublicDocs/EN/CMID/FRI-v1-Packaged%20Product%20-%20Table%20of%20Distributed%20Packages.pdf
# photo year is between 1996 and 2000 
friCH <- st_read("inputNV/ChurchillFRI/FMU_120trlake_175caribou_702lacseul/FMU_120trlake_175caribou_702lacseul.shp")

# Fire AFFES #================================================================
# Use 1980 as oldest fires considered since the data Rempel used included
# "pre-1990" disturbances which meant still visible in 1990 and documentation
# says expect them to be visible ~ 10 years

fireAFFES <- st_read("inputNV/Disturbance/FIREDSTB/LIO-2020-10-28/FIRE_DISTURBANCE_AREA.shp") %>% 
  st_transform(crsUse)

fireAFFESCH2020 <- st_filter(fireAFFES, cHill) %>% 
  filter(between(FIRE_YEAR, 1990, 2020))

fireAFFESCH2010 <- st_filter(fireAFFES, cHill) %>% 
  filter(between(FIRE_YEAR, 1980, 2010))

fireAFFESROF2020 <- st_filter(fireAFFES, rof) %>% 
  filter(between(FIRE_YEAR, 1990, 2020))

fireAFFESROF2010 <- st_filter(fireAFFES, rof) %>% 
  filter(between(FIRE_YEAR, 1980, 2010))

fireAFFESBSand2020 <- st_filter(fireAFFES, bSand) %>% 
  filter(between(FIRE_YEAR, 1990, 2020))

fireAFFESBSand2010 <- st_filter(fireAFFES, bSand) %>% 
  filter(between(FIRE_YEAR, 1980, 2010))

st_write(fireAFFESCH2020, paste0(outCHill, "fireAFFES2020.shp"))
st_write(fireAFFESROF2020, paste0(outROF, "fireAFFES2020.shp"))
st_write(fireAFFESBSand2020, paste0(outBSand, "fireAFFES2020.shp"))

st_write(fireAFFESCH2010, paste0(outCHill, "fireAFFES2010.shp"))
st_write(fireAFFESROF2010, paste0(outROF, "fireAFFES2010.shp"))
st_write(fireAFFESBSand2010, paste0(outBSand, "fireAFFES2010.shp"))

rm(fireAFFES, fireAFFESCH2020, fireAFFESROF2020, fireAFFESBSand2020, 
   fireAFFESCH2010, fireAFFESROF2010, fireAFFESBSand2010)


# HRFCCan Disturbance #=========================================================
# Harvest
HRFCCanHr <- raster("inputNV/Disturbance/CA_forest_harvest_mask_year_1985_2015/CA_harvest_mask_1985_2015.tif")

HRFCCanHrON <- raster::crop(HRFCCanHr, 
                         caribouRanges %>% st_transform(st_crs(HRFCCanHr)))

HRFCCanHrYr <- raster("inputNV/Disturbance/CA_forest_harvest_mask_year_1985_2015/CA_harvest_year_1985_2015.tif")

HRFCCanHrYrON <- raster::crop(HRFCCanHrYr, 
                         caribouRanges %>% st_transform(st_crs(HRFCCanHrYr)))

HRFCCanHrON2010 <- HRFCCanHrON == 1 & HRFCCanHrYrON <= 110

HRFCCanHrCH2010 <- raster::crop(HRFCCanHrON2010, cHill %>% st_transform(st_crs(HRFCCanTy)),
                               filename = paste0(outCHill, "harvHRFCCan2010.tif"),
                               datatype = "INT1U")

HRFCCanHrROF2010 <- raster::crop(HRFCCanHrON2010, rof %>% st_transform(st_crs(HRFCCanTy)),
                                filename = paste0(outROF, "harvHRFCCan2010.tif"),
                                datatype = "INT1U")

HRFCCanHrBSand2010 <- raster::crop(HRFCCanHrON2010, bSand %>% st_transform(st_crs(HRFCCanTy)),
                                  filename = paste0(outBSand, "harvHRFCCan2010.tif"),
                                  datatype = "INT1U")
# 30 years 1985-2015
HRFCCanHrON2015 <- HRFCCanHrON 

HRFCCanHrCH2015 <- raster::crop(HRFCCanHrON2015, cHill %>% st_transform(st_crs(HRFCCanTy)),
                               filename = paste0(outCHill, "harvHRFCCan2015.tif"),
                               datatype = "INT1U")

HRFCCanHrROF2015 <- raster::crop(HRFCCanHrON2015, rof %>% st_transform(st_crs(HRFCCanTy)),
                                filename = paste0(outROF, "harvHRFCCan2015.tif"),
                                datatype = "INT1U")

HRFCCanHrBSand2015 <- raster::crop(HRFCCanHrON2015, bSand %>% st_transform(st_crs(HRFCCanTy)),
                                  filename = paste0(outBSand, "harvHRFCCan2015.tif"),
                                  datatype = "INT1U")

beepr::beep()

# fire
HRFCCanFr <- raster("inputNV/Disturbance/CA_forest_wildfire_year_DNBR_Magnitude_1985_2015/CA_forest_wildfire_year_DNBR_Magnitude_1985_2015.tif",
                    band = 1)

HRFCCanFrON <- raster::crop(HRFCCanFr, 
                            caribouRanges %>% st_transform(st_crs(HRFCCanFr)))

HRFCCanFrYr <- raster("inputNV/Disturbance/CA_forest_wildfire_year_DNBR_Magnitude_1985_2015/CA_forest_wildfire_year_DNBR_Magnitude_1985_2015.tif",
                      band = 2)

HRFCCanFrYrON <- raster::crop(HRFCCanFrYr, 
                              caribouRanges %>% st_transform(st_crs(HRFCCanFrYr)))

HRFCCanFrON2010 <- HRFCCanFrON == 1 & HRFCCanFrYrON <= 110

HRFCCanFrCH2010 <- raster::crop(HRFCCanFrON2010, cHill %>% st_transform(st_crs(HRFCCanFr)),
                                filename = paste0(outCHill, "fireHRFCCan2010.tif"),
                                datatype = "INT1U")

HRFCCanFrROF2010 <- raster::crop(HRFCCanFrON2010, rof %>% st_transform(st_crs(HRFCCanFr)),
                                 filename = paste0(outROF, "fireHRFCCan2010.tif"),
                                 datatype = "INT1U")

HRFCCanFrBSand2010 <- raster::crop(HRFCCanFrON2010, bSand %>% st_transform(st_crs(HRFCCanFr)),
                                   filename = paste0(outBSand, "fireHRFCCan2010.tif"),
                                   datatype = "INT1U")
# 30 years 1985-2015
HRFCCanFrON2015 <- HRFCCanFrON 

HRFCCanFrCH2015 <- raster::crop(HRFCCanFrON2015, cHill %>% st_transform(st_crs(HRFCCanFr)),
                                filename = paste0(outCHill, "fireHRFCCan2015.tif"),
                                datatype = "INT1U")

HRFCCanFrROF2015 <- raster::crop(HRFCCanFrON2015, rof %>% st_transform(st_crs(HRFCCanFr)),
                                 filename = paste0(outROF, "fireHRFCCan2015.tif"),
                                 datatype = "INT1U")

HRFCCanFrBSand2015 <- raster::crop(HRFCCanFrON2015, bSand %>% st_transform(st_crs(HRFCCanFr)),
                                   filename = paste0(outBSand, "fireHRFCCan2015.tif"),
                                   datatype = "INT1U")

beepr::beep()
# CanLaD Disturbance #==========================================================

CanLadYr <- raster("inputNV/Disturbance/CanLaD/CanLaD_20151984_latest_YRT2.tif")
CanLadYrON <- raster::crop(CanLadYr, 
                                caribouRanges %>% st_transform(st_crs(CanLadYr)),
                                filename = "inputNV/Disturbance/CanLaD/CanLaDYrON.tif")

CanLadTy <- raster("inputNV/Disturbance/CanLaD/CanLaD_20151984_latest_type.tif")

CanLadTyON <- raster::crop(CanLadTy, 
                           caribouRanges %>% st_transform(st_crs(CanLadTy)),
                           filename = "inputNV/Disturbance/CanLaD/CanLaDTyON.tif")

# fire
CanLadFrON2010 <- CanLadTyON == 1 & (CanLadYrON <= 2010)

CanLadFrCH2010 <- raster::crop(CanLadFrON2010, cHill %>% st_transform(st_crs(CanLadTy)),
                               filename = paste0(outCHill, "fireCanLad2010.tif"),
                               datatype = "INT1U")

CanLadFrROF2010 <- raster::crop(CanLadFrON2010, rof %>% st_transform(st_crs(CanLadTy)),
                               filename = paste0(outROF, "fireCanLad2010.tif"),
                               datatype = "INT1U")

CanLadFrBSand2010 <- raster::crop(CanLadFrON2010, bSand %>% st_transform(st_crs(CanLadTy)),
                                filename = paste0(outBSand, "fireCanLad2010.tif"),
                                datatype = "INT1U")
# 30 years 1984-2015
CanLadFrON2015 <- CanLadTyON == 1 

CanLadFrCH2015 <- raster::crop(CanLadFrON2015, cHill %>% st_transform(st_crs(CanLadTy)),
                               filename = paste0(outCHill, "fireCanLad2015.tif"),
                               datatype = "INT1U")

CanLadFrROF2015 <- raster::crop(CanLadFrON2015, rof %>% st_transform(st_crs(CanLadTy)),
                                filename = paste0(outROF, "fireCanLad2015.tif"),
                                datatype = "INT1U")

CanLadFrBSand2015 <- raster::crop(CanLadFrON2015, bSand %>% st_transform(st_crs(CanLadTy)),
                                  filename = paste0(outBSand, "fireCanLad2015.tif"),
                                  datatype = "INT1U")

# Harvest
CanLadHrON2010 <- CanLadTyON == 2 & (CanLadYrON <= 2010)

CanLadHrCH2010 <- raster::crop(CanLadHrON2010, cHill %>% st_transform(st_crs(CanLadTy)),
                               filename = paste0(outCHill, "harvCanLad2010.tif"),
                               datatype = "INT1U")

CanLadHrROF2010 <- raster::crop(CanLadHrON2010, rof %>% st_transform(st_crs(CanLadTy)),
                                filename = paste0(outROF, "harvCanLad2010.tif"),
                                datatype = "INT1U")

CanLadHrBSand2010 <- raster::crop(CanLadHrON2010, bSand %>% st_transform(st_crs(CanLadTy)),
                                  filename = paste0(outBSand, "harvCanLad2010.tif"),
                                  datatype = "INT1U")
# 30 years 1984-2015
CanLadHrON2015 <- CanLadTyON == 2 

CanLadHrCH2015 <- raster::crop(CanLadHrON2015, cHill %>% st_transform(st_crs(CanLadTy)),
                               filename = paste0(outCHill, "harvCanLad2015.tif"),
                               datatype = "INT1U")

CanLadHrROF2015 <- raster::crop(CanLadHrON2015, rof %>% st_transform(st_crs(CanLadTy)),
                                filename = paste0(outROF, "harvCanLad2015.tif"),
                                datatype = "INT1U")

CanLadHrBSand2015 <- raster::crop(CanLadHrON2015, bSand %>% st_transform(st_crs(CanLadTy)),
                                  filename = paste0(outBSand, "harvCanLad2015.tif"),
                                  datatype = "INT1U")

beepr::beep()


# Convert Churchill data to 50 m res #==========================================

# plc
plcD <- raster(paste0(outCHill, "plc.tif"))

tmplt_rast <- raster::projectExtent(plcD, crsUseR) %>% raster::`res<-`(50)

plcD <- raster::projectRaster(plcD, tmplt_rast, method = "ngb",
                              filename = paste0(outCHill, "plc50.tif"))

# fire AFFES
fireAFFES2020 <- st_read(paste0(outCHill, "fireAFFES2020.shp"))

fireAFFES2020 <- fasterize::fasterize(fireAFFES2020, tmplt_rast, background = 0)

raster::writeRaster(fireAFFES2020, paste0(outCHill, "fireAFFES2020_50.tif"))

fireAFFES2010 <- st_read(paste0(outCHill, "fireAFFES2010.shp"))

fireAFFES2010 <- fasterize::fasterize(fireAFFES2010, tmplt_rast, background = 0)

raster::writeRaster(fireAFFES2010, paste0(outCHill, "fireAFFES2010_50.tif"))

# HRFCCan
harvHRFCCan2015 <- raster(paste0(outCHill, "harvHRFCCan2015.tif"))

harvHRFCCan2015 <- raster::projectRaster(harvHRFCCan2015, tmplt_rast, 
                                         method = "ngb",
                                     filename = paste0(outCHill,
                                                       "harvHRFCCan2015_50.tif"))

harvHRFCCan2010 <- raster(paste0(outCHill, "harvHRFCCan2010.tif"))

harvHRFCCan2010 <- raster::projectRaster(harvHRFCCan2010, tmplt_rast, 
                                         method = "ngb",
                                     filename = paste0(outCHill,
                                                       "harvHRFCCan2010_50.tif"))

fireHRFCCan2015 <- raster(paste0(outCHill, "fireHRFCCan2015.tif"))

fireHRFCCan2015 <- raster::projectRaster(fireHRFCCan2015, tmplt_rast, 
                                         method = "ngb",
                                         filename = paste0(outCHill,
                                                           "fireHRFCCan2015_50.tif"))

fireHRFCCan2010 <- raster(paste0(outCHill, "fireHRFCCan2010.tif"))

fireHRFCCan2010 <- raster::projectRaster(fireHRFCCan2010, tmplt_rast, 
                                         method = "ngb",
                                         filename = paste0(outCHill,
                                                           "fireHRFCCan2010_50.tif"))

# CanLAD
harvCanLad2015 <- raster(paste0(outCHill, "harvCanLad2015.tif"))

harvCanLad2015 <- raster::projectRaster(harvCanLad2015, tmplt_rast, 
                                         method = "ngb",
                                         filename = paste0(outCHill,
                                                           "harvCanLad2015_50.tif"))

harvCanLad2010 <- raster(paste0(outCHill, "harvCanLad2010.tif"))

harvCanLad2010 <- raster::projectRaster(harvCanLad2010, tmplt_rast, 
                                         method = "ngb",
                                         filename = paste0(outCHill,
                                                           "harvCanLad2010_50.tif"))

fireCanLad2015 <- raster(paste0(outCHill, "fireCanLad2015.tif"))

fireCanLad2015 <- raster::projectRaster(fireCanLad2015, tmplt_rast, 
                                         method = "ngb",
                                         filename = paste0(outCHill,
                                                           "fireCanLad2015_50.tif"))

fireCanLad2010 <- raster(paste0(outCHill, "fireCanLad2010.tif"))

fireCanLad2010 <- raster::projectRaster(fireCanLad2010, tmplt_rast, 
                                         method = "ngb",
                                         filename = paste0(outCHill,
                                                           "fireCanLad2010_50.tif"))


friD <- st_read("inputNV/ChurchillFRI/FMU_120trlake_175caribou_702lacseul/FMU_120trlake_175caribou_702lacseul.shp")
ageD <- fasterize::fasterize(friD, tmplt_rast, field = "AGE")

friLUD <- tibble(PLANFU = friD$PLANFU %>% levels(),
                 code = 1:length(PLANFU)) %>%
  left_join(read.csv("inputNV/ChurchillFRI/OLTBorealCaribouFURegionalForestUnitMapping.csv"),
            by = c(PLANFU = "FU")) %>%
  select(code, Regional.Forest.Unit) %>%
  mutate(Regional.Forest.Unit = toupper(Regional.Forest.Unit))
write.csv(friLUD,
          paste0(outCHill, "friLookUp.csv"), row.names = FALSE)

friD <- fasterize::fasterize(friD, tmplt_rast, field = "PLANFU")
raster::writeRaster(friD, paste0(outCHill, "fri50.tif"))
raster::writeRaster(ageD, paste0(outCHill, "age50.tif"))
