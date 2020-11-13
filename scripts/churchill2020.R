# Current Situation in Churchill

devtools::load_all(".")

# load data
pth_base <- "inputNV/ChurchillData/"

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

coefTableHR %>% filter(Range == "Churchill") %>% 
  pivot_wider(names_from = Season, values_from = Coefficient)
