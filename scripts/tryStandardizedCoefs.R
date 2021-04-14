# Compare results using coefTableStd to results in Hornseth and Rempel 2016

library(caribouMetrics)
library(sf)

caribouRanges <- sf::st_read("../ChurchillAnalysis/inputNV/caribouRanges/Caribou_Range_Boundary.shp", 
                             quiet = TRUE)

esker <- sf::st_read("../ChurchillAnalysis/inputNV/ROFData/esker.shp")
rail <- sf::st_read("../ChurchillAnalysis/inputNV/ROFData/rail.shp")
road_2010 <- sf::st_read("../ChurchillAnalysis/inputNV/ROFData/road_ORNMNRFROF2010.shp")
util_2010 <- sf::st_read("../ChurchillAnalysis/inputNV/ROFData/util2010.shp")


#Bring in the tifs
plcROF <- raster::raster("../ChurchillAnalysis/inputNV/ROFData/plc50.tif") 
fire_AFFES_2010 <- raster::raster("../ChurchillAnalysis/inputNV/ROFData/fireAFFES2010_50.tif")
harv_MNRF_2010 <- raster::raster("../ChurchillAnalysis/inputNV/ROFData/harvMNRF2010_50.tif")

nip <- caribouRanges %>%
  filter(RANGE_NAME %in% c("Nipigon")) %>%
  sf::st_transform(sf::st_crs(plcROF))

missisa <- caribouRanges %>%
  filter(RANGE_NAME %in% c("Missisa")) %>%
  sf::st_transform(sf::st_crs(plcROF)) 

jamesbay <- caribouRanges %>%
  filter(RANGE_NAME %in% c("James Bay")) %>%
  sf::st_transform(sf::st_crs(plcROF)) 

pagwachuan <- caribouRanges %>%
  filter(RANGE_NAME %in% c("Pagwachuan")) %>%
  sf::st_transform(sf::st_crs(plcROF)) 

plcROF <- plcROF %>% reclassPLC()

carHab_Pag <- caribouHabitat(
  landCover = plcROF, 
  esker = esker, 
  natDist = fire_AFFES_2010,  
  anthroDist = harv_MNRF_2010,
  linFeat = road_2010, 
  projectPoly = pagwachuan,
  caribouRange = "Pagwachuan",
  coefTable = coefTableStd,
  doScale = TRUE
)
plot(carHab_Pag)
plot(carHab_Pag, raster.style = "sd", raster.n = 6)
sd_Pag <- raster::cellStats(carHab_Pag@habitatUse, sd)
mean_Pag <- raster::cellStats(carHab_Pag@habitatUse, mean)
# compare to Fig 4 of Hornseth and Rempel 2016 but note that their Summer image
# appears to be a copy of spring
tmap::qtm(raster::scale(carHab_Pag@habitatUse), raster.midpoint = NA,
          raster.breaks = c(-3, -1.5, -0.5, 0.5, 1.5, 2.5), raster.palette = "seq")



