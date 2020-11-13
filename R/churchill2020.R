# Current Situation in Churchill

devtools::load_all(".")

# load data
pth_base <- "inputNV/ChurchillData/"

eskerD <- st_read(paste0(pth_base, "esker.shp"))

# Finish set up only do once convert data to 50 m res
# use the crs of esker data MNR_Lambert_Conformal_Conic 
# crsUse <- st_crs(eskerD) 
# crsUseR <- crs(eskerD)
# 
# plcD <- raster(paste0(pth_base, "plc.tif"))  
# 
# tmplt_rast <- raster::projectExtent(plcD, crsUseR) %>% raster::`res<-`(50)
# 
# plcD <- raster::projectRaster(plcD, tmplt_rast, method = "ngb", 
#                               filename = paste0(pth_base, "plc50.tif"))
# 
# natDistD <- st_read(paste0(pth_base, "fireAFFES2020.shp"))
# 
# natDistD <- fasterize::fasterize(natDistD, tmplt_rast, background = 0)
# 
# raster::writeRaster(natDistD, paste0(pth_base, "fireAFFES2020_50.tif"))
# 
# # change to MNRF Harvest database once aquired
# anthroDistD <- raster(paste0(pth_base, "harvHRFCCan2015.tif")) 
# 
# anthroDistD <- raster::projectRaster(anthroDistD, tmplt_rast, method = "ngb", 
#                                      filename = paste0(pth_base, 
#                                                        "harvHRFCCan2015_50.tif"))
# 
# 
# friD <- st_read("inputNV/ChurchillFRI/FMU_120trlake_175caribou_702lacseul/FMU_120trlake_175caribou_702lacseul.shp")
# ageD <- fasterize::fasterize(friD, tmplt_rast, field = "AGE") 
# 
# friLUD <- tibble(PLANFU = friD$PLANFU %>% levels(), 
#                  code = 1:length(PLANFU)) %>% 
#   left_join(read.csv("inputNV/ChurchillFRI/OLTBorealCaribouFURegionalForestUnitMapping.csv"),
#             by = c(PLANFU = "FU")) %>% 
#   select(code, Regional.Forest.Unit) %>% 
#   mutate(Regional.Forest.Unit = toupper(Regional.Forest.Unit))
# write.csv(friLUD,
#           paste0(pth_base, "friLookUp.csv"), row.names = FALSE)
# 
# friD <- fasterize::fasterize(friD, tmplt_rast, field = "PLANFU") 
# raster::writeRaster(friD, paste0(pth_base, "fri50.tif"))
# raster::writeRaster(ageD, paste0(pth_base, "age50.tif"))

eskerD <- st_read(paste0(pth_base, "esker.shp"))
plcD <- raster(paste0(pth_base, "plc50.tif"))
roadD <- st_read(paste0(pth_base, "road_ORNMNRFCH2020.shp"))
railD <- st_read(paste0(pth_base, "rail.shp"))
utilD <- st_read(paste0(pth_base, "util2020.shp"))
natDistD <- raster(paste0(pth_base, "fireAFFES2020_50.tif"))
anthroDistD <- raster(paste0(pth_base, "harvHRFCCan2015_50.tif"))
friD <- raster(paste0(pth_base, "fri50.tif"))
ageD <- raster(paste0(pth_base, "age50.tif"))
friLUD <- read.csv(paste0(pth_base, "friLookUp.csv"))

projectPolyD <- st_read("./inputNV/caribouRanges/Caribou_Range_Boundary.shp", 
                        quiet = TRUE) %>% 
  filter(RANGE_NAME == "Churchill") %>% 
  st_transform(crs = st_crs(plcD))

caribouResults <- caribouHabitat(
  plc = plcD, 
  esker = eskerD, 
  fri = friD, 
  age = ageD, 
  natDist = natDistD, 
  anthroDist = anthroDistD,
  harv = anthroDistD,
  # linFeat = linFeatDen,
  linFeat = list(roads = roadD, rail = railD,
                 utilities = utilD),
  projectPoly = projectPolyD, 
  friLU = friLUD, 
  caribouRange = "Churchill", 
  padProjPoly = TRUE, 
  eskerSave = paste0(pth_base, "eskerDen.tif"),
  linFeatSave = paste0(pth_base, "linFeatDen.tif")
) 
beepr::beep()

caribouResults@fri
caribouResults@plc

caribouResults@processedData

plot(caribouResults)

tmap::qtm(caribouResults@processedData[[c(1:6, 8:12)]] %>% 
            raster::mask(caribouResults@projectPoly),
          facets.free.scales = TRUE, 
          layout.legend.show = FALSE)

coefTableHR %>% filter(Range == "Churchill") %>% arrange(Variable)
