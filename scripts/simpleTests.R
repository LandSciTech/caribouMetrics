require(raster)
require()
cWidth=300

#############
# Simple polygon dataset test
#With buffers, mismatch between theory and calculations is <=0.1%
testData =makeTestPolygons(cWidth)
class(testData$cLS)
#raster::plot(testData$cLS@maps$ranges)
#raster::plot(testData$cLS@maps$fires[[1]],add=T,col="red")
#raster::plot(testData$cLS@maps$anthroPoly[[1]],add=T,col="blue")
#raster::plot(testData$cLS@maps$linear[[1]],add=T,col="black")

identical(round(testData$answersBuffer,1),
  round(subset(bCaribouMetaAnalysisPredictors(testData$cLS,width=cWidth),select=c(fire,anthro,totalDist)),1))

identical(round(testData$answersNoBuffer,6),round(subset(bCaribouMetaAnalysisPredictors(testData$cLS,width=0),select=c(fire,anthro,totalDist)),6))

identical(round(testData$answersBufferHalf,1),
          round(subset(bCaribouMetaAnalysisPredictors(testData$cLS,width=cWidth*0.5),select=c(fire,anthro,totalDist)),1))

##############
# Test rasters
#accuracy will depend on discretization. 1.35% for this example.
testData = makeTestPolygons(cWidth,makeRaster=c("ranges","anthroPoly","fires"))

#note this is slow
identical(round(testData$answersBuffer,2),
          round(subset(bCaribouMetaAnalysisPredictors(testData$cLS,width=cWidth),select=c(fire,anthro,totalDist)),2))

identical(round(testData$answersNoBuffer,2),round(subset(bCaribouMetaAnalysisPredictors(testData$cLS,width=0),select=c(fire,anthro,totalDist)),2))

checkHalf = bCaribouMetaAnalysisPredictors(testData$cLS,width=cWidth*0.5)
checkHalf= subset(checkHalf, select=c(fire,anthro,totalDist))
max(testData$answersBufferHalf - checkHalf)<1.35
rasterDat = testData$cLS

#########
#Test raster conversion routine
testData = makeTestPolygons(cWidth,makeRaster=c("ranges"),forceRaster=T)
identical(bCaribouMetaAnalysisPredictors(rasterDat,width=0),bCaribouMetaAnalysisPredictors(testData$cLS,width=0))

################
#TO DO test raster/polygon combinations.
#TO DO build formal test suite.
#TO DO start to construct vignette
