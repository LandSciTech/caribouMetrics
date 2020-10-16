# Explore the effect of disturbance within model
library(tmap)
library(ggplot2)
library(purrr)
devtools::load_all(".")

# load data
sSimData <- getSyncSimData(ssimLib = "./inputNV/Libraries/ChurchillBC",
                           ssimProject = "2008/09 Plans",
                           rfuLUPath = "inputTables/OLTBorealCaribouFURegionalForestUnitMapping.csv",
                           scnNum = 6)

plcD <- raster("./inputNV/intermediaryData/plc_aligned.tif")
eskerD <- raster("inputNV/intermediaryData/eskerRastDen.tif")
natDistD <- raster("./inputNV/intermediaryData/natDistOFRI.tif")
anthroDistD <- raster("./inputNV/intermediaryData/dummyAnthroDist.tif")
harvD <- raster("inputNV/intermediaryData/harvPost2000OFRI.tif")
linFeatD <- raster("./inputNV/intermediaryData/linFeatRastDen.tif")
projectPolyD <- st_read("./inputNV/caribouRanges/Caribou_Range_Boundary.shp", 
                        quiet = TRUE) %>% 
  filter(RANGE_NAME == "Churchill")

# Normal version
resWithDist <- caribouHabitat(plc = plcD, esker = eskerD, fri = sSimData$mySC, 
                              age = sSimData$myAge, natDist = natDistD, 
                              anthroDist = anthroDistD,
                              harv = harvD,
                              linFeat = linFeatD, projectPoly = projectPolyD, 
                              friLU = sSimData$friLU, 
                              caribouRange = "Churchill")

# Change the plc data so that it is all natural no disturbance
natDistD0 <- raster::setValues(natDistD, 0)
anthroDistD0 <- raster::setValues(anthroDistD, 0)
harvD0 <- raster::setValues(harvD, 0)

plcD %>% raster::values() %>% table() %>% as.data.frame() %>% 
  mutate(Code = as.integer(as.character(`.`))) %>% 
  right_join(plcToResType, by = c(Code = "PLCCode"))

plcToResType <- plcToResType %>% 
  mutate(ResourceType = case_when(PLCCode == 7 ~ "other",
                                  PLCCode == 8 ~ "other",
                                  PLCCode == 9 ~ "other", 
                                  TRUE ~ ResourceType))

# No disturbance version
resNoDist <- caribouHabitat(plc = plcD, esker = eskerD, fri = sSimData$mySC, 
                              age = sSimData$myAge, natDist = natDistD0, 
                              anthroDist = anthroDistD0,
                              harv = harvD0,
                              linFeat = linFeatD, projectPoly = projectPolyD, 
                              friLU = sSimData$friLU, 
                              caribouRange = "Churchill")

resNoDist@habitatUse %>% plot() 
resWithDist@habitatUse %>% plot()

resNoDist@plc %>% tmapCatRast(lbls = raster::levels(resNoDist@plc)[[1]]$resType,
                              ttle = "NoDist")

resWithDist@plc %>% tmapCatRast(lbls = raster::levels(resWithDist@plc)[[1]]$resType,
                              ttle = "WithDist")

resNoDist@processedData$CON %>% plot(main = "NoDist CON")
resWithDist@processedData$CON %>% plot(main = "WithDist CON")

binNoDist <- calcBinaryUse(resNoDist, caribouRange = "Churchill") 
binNoDist %>% plot(main = "NoDist Cat2")

binWithDist <- calcBinaryUse(resWithDist, caribouRange = "Churchill") 
binWithDist %>% plot(main = "WithDist Cat2")

binNoDist %>% raster::cellStats(sum) * raster::res(binNoDist)[1]^2 / 10000 - 
  binWithDist %>% raster::cellStats(sum) * raster::res(binWithDist)[1]^2 / 10000

# The difference between no disturbance and with disturbance is 434400 ha or
# 11.7 % of the Cat2 area with no disturbance General patterns are similar and
# are mostly driven by landcover as opposed to disturbance

coefGrph <- coefTableHR %>% 
  ggplot(aes(Variable, Coefficient, fill = Season))+
  geom_col(position = "dodge")+
  coord_flip()+
  facet_wrap(~Range)

coefTableHR %>% filter(Range %in% c("Missisa", "James Bay", "Ozhiski")) %>% 
  ggplot(aes(Variable, Coefficient, fill = Season))+
  geom_col(position = "dodge")+
  coord_flip()+
  facet_wrap(~Range)

# NOTE: that Missisa (Where the ROF is) has by far the strongest effect of linear
# features this may be related to the fact it has the lowest amount of linear
# features

