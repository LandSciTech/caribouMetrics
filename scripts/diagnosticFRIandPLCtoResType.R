# Diagnostic on the conversion from FRI and PLC to Resource types 
library(tmap)
library(ggplot2)
devtools::load_all(".")


sSimData <- getSyncSimData(ssimLib = "./inputNV/Libraries/ChurchillBC",
                           ssimProject = "2008/09 Plans",
                           rfuLUPath = "inputTables/OLTBorealCaribouFURegionalForestUnitMapping.csv",
                           scnNum = 6)


# data as rasters and hexgrid provided
plcD <- raster("./inputNV/intermediaryData/plc_aligned.tif")
#eskerD <- st_read("inputNV/HornsethRempelInfo/Eskers_Ontario/Eskers_Ontario.shp")
eskerD <- raster("inputNV/intermediaryData/eskerRastDen.tif")
natDistD <- raster("./inputNV/intermediaryData/dummyNatDist.tif")
anthroDistD <- raster("./inputNV/intermediaryData/dummyAnthroDist.tif")
linFeatD <- raster("./inputNV/intermediaryData/linFeatRastDen.tif")
#linFeatD <- st_read("inputNV/Roads/TroutLakeRoads_ornsegad.shp")
projectPolyD <- st_read("./inputNV/intermediaryData/projectPoly.shp")
hexgridD <- st_read("./inputNV/intermediaryData/testPoly.shp")

caribouResults <- caribouHabitat(plc = plcD, esker = eskerD, fri = sSimData$mySC, 
                                 age = sSimData$myAge, natDist = natDistD, 
                                 anthroDist = anthroDistD,
                                 linFeat = linFeatD, projectPoly = projectPolyD, 
                                 friLU = sSimData$friLU, 
                                 caribouRange = "Churchill", 
                                 hexgrid = hexgridD)

# pick a subset that shows area that overlaps boundary of area covered by fri
# plot(caribouResults@fri)
# ext <- drawExtent()
ext <- raster::extent(c(xmin = 334689.6, xmax = 382137.2,
                        ymin = 12676724, ymax = 12715753))
 
croppedData <- lst(caribouResults@plc, caribouResults@fri, caribouResults@age, 
                   plcD, sSimData$mySC) %>% 
  purrr::map(~raster::crop(.x, ext))

# get resType names
nms <- croppedData$`caribouResults@plc` %>% raster::levels() %>% .[[1]] %>% 
  pull(resType) %>% as.character()

# RAT messes up plotting
croppedData$`caribouResults@plc` <- croppedData$`caribouResults@fri` %>% 
  raster::deratify(att = 1, drop = TRUE, complete = TRUE) 

tmap_mode("plot")

plotResType <- tmapCatRast

plcMap <- plotResType(croppedData$`caribouResults@plc`, nms, "PLC ResType")

friMap <- plotResType(croppedData$`caribouResults@fri`, nms, "FRI ResType")

# histogram comparing PLC and FRI resource types
friPlcHist <- raster::freq(caribouResults@fri) %>% as.data.frame() %>% 
  filter(!is.na(value)) %>% 
  mutate(resType = nms[.$value], 
         dataSource = "FRI") %>%
  bind_rows(raster::freq(caribouResults@plc) %>% as.data.frame() %>% 
              filter(!is.na(value)) %>% 
              mutate(resType = nms, dataSource = "PLC")) %>% 
  ggplot(aes(resType, count, fill = dataSource))+
  geom_col(position = "dodge")+
  expand_limits(y = 4e+05)+
  labs(x = "Resource Type")+
  theme_classic()


# adjust rfuToResType so that OcLow and SbLow are labeled ST and PjMx1 and SbMx1
# are labeled MIX
rfuToResType <- mutate(rfuToResType, 
                       ResourceType = case_when(
                         RegionalForestUnit == "SbLow" ~ "ST",
                         RegionalForestUnit == "OcLow" ~ "ST",
                         RegionalForestUnit == "PjMx1" ~ "MIX",
                         RegionalForestUnit == "SbMx1" ~ "MIX",
                         TRUE ~ ResourceType)
                       ) 

caribouResults2 <- caribouHabitat(plc = plcD, esker = eskerD, fri = sSimData$mySC, 
                                 age = sSimData$myAge, natDist = natDistD, 
                                 linFeat = linFeatD, projectPoly = projectPolyD, 
                                 friLU = sSimData$friLU, 
                                 caribouRange = "Churchill", 
                                 hexgrid = hexgridD)
beepr::beep()
# histogram comparing PLC and FRI resource types
friPlcHist2 <- raster::freq(caribouResults2@fri) %>% as.data.frame() %>% 
  filter(!is.na(value)) %>% 
  mutate(resType = nms, 
         dataSource = "FRI") %>%
  bind_rows(raster::freq(caribouResults2@plc) %>% as.data.frame() %>% 
              filter(!is.na(value)) %>% 
              mutate(resType = nms, dataSource = "PLC")) %>% 
  ggplot(aes(resType, count, fill = dataSource))+
  geom_col(position = "dodge")+
  expand_limits(y = 4e+05)+
  labs(x = "Resource Type")+
  theme_classic()

# adding SbLow and OcLow to ST and adding PjMx1 and SbMx1 to MIX improved
# correspondence
friPlcHist2

# NOTE: so far no change made to stored version. Rob Rempel said that the
# crosswalk was designed looking at the whole area of the PLC so changes to
# match one FMU might lead to a poorer match for other areas.