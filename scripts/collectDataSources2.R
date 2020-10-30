# create data sets for Churchill and ROF

devtools::load_all(".")

outCHill <- "inputNV/ChurchillData/"
outROF <- "inputNV/ROFData/"

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

st_write(eskerCH, paste0(outCHill, "esker.shp"))
st_write(eskerROF, paste0(outROF, "esker.shp"))

rm(eskerD, eskerCH, eskerROF)

# road #========================================================================

# Need 2010 and 2020 versions use year constructed in MNRF data
road_ORN <- st_read("inputNV/Roads/ORNSEGAD_20201028/Non_Sensitive.gdb")

road_MNRF <- st_read("inputNV/Roads/MNRRDSEG_20201028/Non_Sensitive.gdb", 
                     layer = "MNRF_ROAD_SEGMENT")

road_ORNCH <- st_filter(road_ORN, cHill %>% st_transform(st_crs(road_ORN))) %>% 
  st_transform(crsUse)

road_MNRFCH <- st_filter(road_MNRF, cHill %>% st_transform(st_crs(road_MNRF))) %>% 
  st_transform(crsUse)

st_write(road_MNRFCH %>% st_cast("LINESTRING"),
         paste0(outCHill, "temproad_MNRFCH.gpkg"), append = FALSE)
st_write(road_ORNCH %>% st_cast("LINESTRING"), 
         paste0(outCHill, "temproad_ORNCH.gpkg"), append = FALSE)

# Used QGIS to buffer MNRF by 100 m and then difference that with ORN. So ORN -
# MNRF and saved as road_ORNMNRF.gpkg

road_ORNMNRFCH <- st_read(paste0(outCHill, "road_ORNMNRFCH.gpkg")) %>% 
  st_set_crs(crsUse) %>% select(YEAR_CONSTRUCTED)

road_ORNMNRFCH2010 <- road_ORNMNRFCH %>% 
  filter(YEAR_CONSTRUCTED <= 2010 | is.na(YEAR_CONSTRUCTED))

plot(road_ORNMNRFCH %>% st_geometry(), col = "red")
plot(road_ORNMNRFCH2010 %>% st_geometry(), add = TRUE)
plot(road_ORNMNRFCH %>% filter(is.na(YEAR_CONSTRUCTED)) %>% st_geometry(), 
     col = "green", add = TRUE)
plot(road_ORNMNRFCH %>% filter(YEAR_CONSTRUCTED == 9999) %>% st_geometry(), 
     col = "blue", add = TRUE)

st_write(road_ORNMNRFCH, paste0(outCHill, "road_ORNMNRFCH2020.shp"))
st_write(road_ORNMNRFCH2010, paste0(outCHill, "road_ORNMNRFCH2010.shp"))
