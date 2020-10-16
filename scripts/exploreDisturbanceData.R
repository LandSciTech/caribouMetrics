# Explore Data sources
devtools::load_all(".")
library(purrr)
library(tmap)

# Explore OFRI disturbance data #===============================================
# Load series of points
basePth <- "./inputNV/Disturbance/OFRIDisturbanceHistory/ChurchillCaribouBorealDisturbance"

pths <- list.files(basePth, pattern = ".shp$", full.names = TRUE)

# load all the points into one df
OFRIDistHist <- map(pths, st_read) %>% rbind_list()

OFRIDist1972 <- pths[1] %>% st_read()
# Points are only present if Disturbed = 1. The same point may be present in
# multiple years, identified by pointid

# Look at timeseries for one point
OFRIDistHist %>% filter(pointid == OFRIDistHist$pointid[1000001]) %>% View()

# Look at all timeseries
OFRIDistHist %>% arrange(pointid, Disturbanc) %>% View()

# Natural disturbance if it was natDist and is >35 then it will be updated Don't
# really know how long before PLC in 2000 to consider? The PLC manual says that
# Burns are defined as burned in the last 10 years
natDistOFRI <- OFRIDistHist %>% 
  filter(FireConfid > 0, Disturbanc > 1990) %>% 
  group_by(pointid) %>% 
  summarise(Disturbanc = first(Disturbanc), 
            geometry = first(geometry)) %>% 
  st_as_sf() %>% 
  st_set_crs(st_crs(OFRIDist1972)) %>% 
  st_transform(st_crs(plcD))

# Harvest Disturbance is used to determine whether the area was harvested after
# the PLC in 2000 was made and so should be updated with FRI
harvPost2000OFRI <- OFRIDistHist %>% 
  filter(HarvestCon > 0, Disturbanc > 2000) %>% 
  group_by(pointid) %>% 
  summarise(Disturbanc = first(Disturbanc), 
            geometry = first(geometry)) %>% 
  st_as_sf() %>% 
  st_set_crs(st_crs(OFRIDist1972)) %>% 
  st_transform(st_crs(plcD))

  
# make it into a raster with 1 for disturbed and 0 otherwise
natDistOFRIRas <- natDistOFRI %>% filter(!st_is_empty(geometry)) %>%
  mutate(geometry = st_cast(geometry, "POINT")) %>% 
  raster::rasterize(plcD, field = 1, background = 0)

raster::writeRaster(natDistOFRIRas, 
                    "inputNV/intermediaryData/natDistOFRI.tif")

harvPost2000OFRIRas <- harvPost2000OFRI %>% filter(!st_is_empty(geometry)) %>%
  mutate(geometry = st_cast(geometry, "POINT")) %>% 
  raster::rasterize(plcD, field = 1, background = 0)

raster::writeRaster(harvPost2000OFRIRas, 
                    "inputNV/intermediaryData/harvPost2000OFRI.tif")
# Visualize rasterizing results
ext <- raster::drawExtent()

natDistOFRIRasCrop <- raster::crop(natDistOFRIRas, ext)
natDistOFRICrop <- st_crop(natDistOFRI, ext)

harvOFRIRasCrop <- raster::crop(harvPost2000OFRIRas, ext)
harvOFRICrop <- st_crop(harvPost2000OFRI, ext)

plot(harvOFRIRasCrop)
plot(st_geometry(harvOFRICrop), add = TRUE)

# Explore Provincial Satellite Derived disturbance data #=======================
# This is the one referenced in the Hornseth Rempel report

gdbPth <- "inputNV/Disturbance/ProvSatelliteDerivedDisturbanceMapping/Provincial_Satellite_Derived_Disturbance_Mapping-FarNorth_BorealShield.gdb"
st_layers(gdbPth)
provSatDist15 <- st_read(gdbPth, layer = "DISTURBANCE_Z15")
provSatDist16 <- st_read(gdbPth, layer = "DISTURBANCE_Z16")
provSatDist15 <- st_make_valid(provSatDist15)

tmap_mode("view")
qtm(st_bbox(provSatDist16) %>% st_as_sfc())

qtm(st_bbox(provSatDist15) %>% st_as_sfc())+ qtm(projectPolyD %>% st_geometry())

provSatNatDist <- provSatDist %>% st_transform(st_crs(projectPolyD)) %>% 
  st_crop(projectPolyD)

beepr::beep()

# Cuts off half way through our range! Does not seem to be available for
# southern part of Ontario?


# Explore wooded area #=========================================================
wooded <- st_read("inputNV/WOODAREA/LIO-2019-09-30/WOODED_AREA.shp")

projectPolyD <- st_read("./inputNV/caribouRanges/Caribou_Range_Boundary.shp") %>% 
  filter(RANGE_NAME == "Churchill")

woodedChill <- st_crop(wooded %>% st_transform(st_crs(projectPolyD)), 
                       projectPolyD)

tmap_mode("view")
tm_shape(woodedChill)+
  tm_fill()

# Doesn't seem very helpful

# Explore original FRI data for Churchill range #===============================
friRaw <- st_read("inputNV/FMU_120trlake_175caribou_702lacseul/FMU_120trlake_175caribou_702lacseul.shp")


friRaw %>% st_drop_geometry() %>%  group_by(SFU) %>% summarise(nSFU = n()) %>% View()

# present seems to be 2012

# Is there a field that represents Forest?

# From Josie: interpretation of the FRI fields is encoded in
# \lstools\R\STSimONFRIInitialConditions.R. For the Churchill model
# non-productive forest is defined using a combination of DEVSTAGE and OWNER.
# non-productive default is DEVSTAGE == "LOWNAT" and OWNER == 1

# Need to know stands that have been harvested since PLC was made ie 2000
harvPost2000 <- friRaw %>%
  filter(YRHARV >= 2000, stringr::str_detect(DEPTYPE, "Harvest"))

harvPost2000Rast <- fasterize::fasterize(harvPost2000, sSimData$mySC)

raster::writeRaster(harvPost2000Rast, 
                    "inputNV/intermediaryData/harvPost2000.tif", overwrite = TRUE)

# Use YRFIRE to find Natural disturbance at time of PLC
friRaw %>% st_drop_geometry() %>% group_by(DEPTYPE) %>% summarise(N = n())

friRaw %>% filter(stringr::str_detect(DEPTYPE, "Fire"), YRFIRE < 2000) %>% View()

# Fire.f and Fire.p are only used for recent fires
firePre2000 <- friRaw %>% filter(YRFIRE < 2000, DEPTYPE == "Fire")

firePre2000Rast <- fasterize::fasterize(firePre2000, sSimData$mySC)

raster::writeRaster(firePre2000Rast, 
                    "inputNV/intermediaryData/firePre2000.tif", overwrite = TRUE)



