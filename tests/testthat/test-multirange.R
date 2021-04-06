
# path to use interactively
# pthBase <- "tests/testthat/data/"
pthBase <- "data/"


landCoverD = raster(paste0(pthBase, "plc", ".tif")) %>% 
  reclassPLC()
eskerDras = raster(paste0(pthBase, "eskerTif", ".tif"))
eskerDshp = st_read(paste0(pthBase, "esker", ".shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
friLUD = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE)
updatedLCD = raster(paste0(pthBase, "fri", ".tif")) %>% 
  reclassFRI(friLUD)
ageD = raster(paste0(pthBase, "age", ".tif"))
natDistD = raster(paste0(pthBase, "natDist", ".tif"))
anthroDistD = raster(paste0(pthBase, "anthroDist", ".tif"))
harvD = raster(paste0(pthBase, "harv", ".tif"))
linFeatDras = raster(paste0(pthBase, "linFeatTif", ".tif"))
projectPolyD = st_read(paste0(pthBase, "projectPoly", ".shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDshp = st_read(paste0(pthBase, "linFeat", ".shp"), quiet = TRUE) %>% 
  st_set_agr("constant")

data_esktif_linFtif <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras, updatedLC = updatedLCD, 
  age = ageD, natDist = natDistD,
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDras, projectPoly = projectPolyD,
  caribouRange = "Churchill", 
  winArea = 500
)

projPolyPts <- projectPolyD$geometry %>% st_cast("POINT") %>% st_coordinates()
twoRange <- st_sfc(st_polygon(list(projPolyPts[c(1, 2, 3 ,1),])), 
                   st_polygon(list(projPolyPts[c(1, 4, 3 ,1),]))) %>%
  st_as_sf() %>% 
  mutate(ID = 1:n(), 
         Range = c("Missisa", "Nipigon")) %>% 
  st_set_crs(st_crs(projectPolyD))

# supply polygon with multiple ranges

# same coefficients as range
resTwoRange <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras, updatedLC = updatedLCD, 
  age = ageD, natDist = natDistD,
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDras, projectPoly = twoRange,
  caribouRange = data.frame(Range = c("Missisa", "Nipigon"), 
                            coefRange = c("Missisa", "Nipigon"), 
                            stringsAsFactors = FALSE), 
  winArea = 500
)

# different coefficients as range
resTwoRangeDif <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras, updatedLC = updatedLCD, 
  age = ageD, natDist = natDistD,
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDras, projectPoly = twoRange,
  caribouRange = data.frame(Range = c("Missisa", "Nipigon"), 
                            coefRange = c("Nipigon", "Missisa"), 
                            stringsAsFactors = FALSE), 
  winArea = 500
)

# different winArea need to supply dif coef table to work on small example
coefTableSmall <- coefTableHR %>% mutate(WinArea = WinArea/10)

resTwoRangeDifWin <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras, updatedLC = updatedLCD, 
  age = ageD, natDist = natDistD,
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDras, projectPoly = twoRange,
  caribouRange = data.frame(Range = c("Missisa", "Nipigon"), 
                            coefRange = c("Nipigon", "Missisa"), 
                            stringsAsFactors = FALSE), 
  coefTable = coefTableSmall
)

resRange1 <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras, updatedLC = updatedLCD, 
  age = ageD, natDist = natDistD,
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDras, projectPoly = twoRange %>% slice(1),
  caribouRange = data.frame(Range = c("Missisa"), 
                            coefRange = c("Nipigon"), 
                            stringsAsFactors = FALSE), 
  coefTable = coefTableSmall
)

resRange2 <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras, updatedLC = updatedLCD, 
  age = ageD, natDist = natDistD,
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDras, projectPoly = twoRange %>% slice(2),
  caribouRange = data.frame(Range = c("Nipigon"), 
                            coefRange = c("Missisa"), 
                            stringsAsFactors = FALSE), 
  coefTable = coefTableSmall
)

# TODO: figure out why these are so different
plot(resTwoRangeDifWin@processedData$CON)
plot(sum(resRange1@processedData$CON, resRange2@processedData$CON, 
         na.rm = TRUE))
