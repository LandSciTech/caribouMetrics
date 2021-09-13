#profile mem use
library(sf)
library(raster)
library(dplyr)
library(tidyr)
library(purrr)

# I compared to old versions by installing the current version of the package
# but not loadsing it and then making changes and sourcing the package functions
# so that the new version is in GlobalEnvironment and old can be accessed with
# caribouMetrics::

pthBase <- "tests/testthat/data/"
plcToResType <- caribouMetrics::plcToResType
resTypeCode <- caribouMetrics::resTypeCode
coefTableHR <- caribouMetrics::coefTableHR
# # load data for tests
# landCoverD = raster(paste0(pthBase, "landCover", ".tif")) %>% 
#   caribouMetrics::reclassPLC()
# eskerDras = raster(paste0(pthBase, "eskerTif", ".tif"))
# eskerDshp = st_read(paste0(pthBase, "esker", ".shp"), quiet = TRUE)
# natDistD = raster(paste0(pthBase, "natDist", ".tif"))
# anthroDistD = raster(paste0(pthBase, "anthroDist", ".tif"))
# linFeatDras = raster(paste0(pthBase, "linFeatTif", ".tif"))
# projectPolyD = st_read(paste0(pthBase, "projectPoly", ".shp"), quiet = TRUE)
# linFeatDshp = st_read(paste0(pthBase, "linFeat", ".shp"), quiet = TRUE)
# 
# memcH <- profmem::profmem({
#   cH <- caribouHabitat(
#     landCover = landCoverD, 
#     esker = eskerDras, 
#     natDist = natDistD, 
#     anthroDist = anthroDistD, 
#     linFeat = linFeatDras, 
#     projectPoly = projectPolyD,
#     caribouRange = "Nipigon",
#     winArea = 500)
# })
# inData <- inputData(
#   landCover = landCoverD, 
#   esker = eskerDras, 
#   natDist = natDistD, 
#   anthroDist = anthroDistD, 
#   linFeat = linFeatDras, 
#   projectPoly = projectPolyD,
#   caribouRange = "Nipigon",
#   winArea = 500,
#   coefTable = coefTableHR)
# 
# procData <- processData(inData)
# 
# # tried using calc
# outCur <- calcRSP(procData@processedData, coefTableHR %>% filter(Range == "Nipigon"))
# outNew <- newCalcRSP(procData@processedData, coefTableHR %>% filter(Range == "Nipigon"))
# # new version is 10 times slower on small dataset
# 
# 
# resourceProp <- procData@processedData
# coefs <- coefTableHR %>% filter(Range == "Nipigon")
# seasons = "all"
# doScale = FALSE
# 
# # apparently better for setting values 
# set1f <- function(x){rep(1, x)}
# z1 <- raster::init(landCoverD, fun=function(x){rep(1, x)}, 
#                    filename = raster::rasterTmpFile())

# changed raster so that uses mask instead of [] <- 
sourceDir <- function(path, trace = TRUE, ...) {
  op <- options(); on.exit(options(op)) # to reset after each 
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
    options(op)
  }
}
sourceDir("R")

# profvis::profvis({
#   dM <- disturbanceMetrics(
#     landCover = landCoverD,  
#     natDist = natDistD, 
#     anthroDist = anthroDistD, 
#     linFeat = linFeatDras, 
#     projectPoly = projectPolyD)
# })
# 
# profvis::profvis({
#   dM2 <- caribouMetrics::disturbanceMetrics(
#     landCover = landCoverD,  
#     natDist = natDistD, 
#     anthroDist = anthroDistD, 
#     linFeat = linFeatDras, 
#     projectPoly = projectPolyD)
# })

# try with larger data
pth_base <- "vignettes/Example_data/"

landCover <- raster(paste0(pth_base, "plc250.tif"))%>%
  reclassPLC()

esker <- st_read(paste0(pth_base, "esker.shp"))

linFeat <- list(roads = st_read(paste0(pth_base, 
                                       "road_ORNMNRFROF2010.shp")),
                rail = st_read(paste0(pth_base, "rail.shp")),
                utilities = st_read(paste0(pth_base, "util2010.shp")))

natDist <- raster(paste0(pth_base, "natDist250.tif"))

projectPoly <- st_read(paste0(pth_base,
                              "caribouRanges/Caribou_Range_Boundary.shp")) %>% 
  filter(RANGE_NAME == "Nipigon") %>% 
  rename(Range = RANGE_NAME)

singlePoly <- projectPoly[1,]

memDMnew3 <- profmem::profmem({
  disturb <- disturbanceMetrics(landCover = landCover,
                                linFeat = linFeat,  
                                natDist = natDist,
                                projectPoly = singlePoly)
})

# memDMold <- profmem::profmem({
#   disturb2 <- caribouMetrics::disturbanceMetrics(landCover = landCover,
#                                                 linFeat = linFeat,  
#                                                 natDist = natDist,
#                                                 projectPoly = singlePoly)
# })
bm <- bench::mark(processData(disturb), caribouMetrics::processData(disturb), 
                  check = FALSE)
bm3 <- bench::mark(processData(disturb), check = FALSE)

# Look at caribouHabitat
tmplt <- raster(disturb@landCover) %>% raster::`res<-`(c(400, 400))
profvis::profvis(applyDist(disturb@landCover, disturb@natDist, disturb@anthroDist, 
                           tmplt))

bm4 <- bench::mark(new = applyDist(disturb@landCover, disturb@natDist, 
                                   disturb@anthroDist, tmplt),
                   old = caribouMetrics::applyDist(disturb@landCover, disturb@natDist, 
                                                   disturb@anthroDist, tmplt))
bm5 <-  bench::mark(new2 = applyDist(disturb@landCover, disturb@natDist, 
                                    disturb@anthroDist, tmplt))

inData <- inputData(
  landCover = landCover,
  esker = esker,
  natDist = natDist,
  linFeat = linFeat,
  projectPoly = singlePoly,
  caribouRange = "Nipigon",
  coefTable = coefTableHR,
  winArea = 5000)

bmprocD <- bench::mark(old = old <- caribouMetrics::processData(inData),
                       new = new <- processData(inData), check = FALSE)

