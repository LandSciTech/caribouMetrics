
devtools::load_all()
pthBase <- "./tests/testthat/data/"

# load example data and classify plc and fri into Resource Types
plcD = raster(paste0(pthBase, "plc", ".tif")) 
plcD[plcD==1]=NA #Assume 1 is water, all other values are irrelevant
natDistD = raster(paste0(pthBase, "natDist", ".tif"))
anthroDistD = raster(paste0(pthBase, "anthroDist", ".tif"))
harvD = raster(paste0(pthBase, "harv", ".tif"))
projectPolyD = st_read(paste0(pthBase, "projectPoly", ".shp"), quiet = TRUE)
linFeatDshp = st_read(paste0(pthBase, "linFeat", ".shp"), quiet = TRUE)
roadsD = st_read(paste0(pthBase, "roads.shp"), quiet = TRUE)
railD = st_read(paste0(pthBase, "rail.shp"), quiet = TRUE)
utilitiesD = st_read(paste0(pthBase, "utilities.shp"), quiet = TRUE)

devtools::load_all()

#Note assumption here is natDist, anthroDist and harv include 40 years of
#cumulative disturbance.
dm <- disturbanceMetrics(
  landCover=plcD,
  natDist = natDistD, 
  anthroDist = anthroDistD, 
  harv = harvD,
  linFeat = linFeatDshp, 
  projectPoly = projectPolyD,
  padFocal = TRUE, # assume data outside area is 0 for all variables
  bufferWidth = 500 # only use in examples, leave as default for correct results
)

dm@disturbanceMetrics
plot(dm@processedData)

#TODO: testing, documentation, speeding up, etc.