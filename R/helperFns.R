#' @include AAAClassDefinitions.R
NULL


getAllYrs<-function(allYrs){
  if(!identical('static',allYrs)){
    allYrs=setdiff(allYrs,'static')
    allYrs = sort(as.numeric(gsub("t","",allYrs,fixed=T)))
  }
  return(allYrs)
}
#' Test function for calculation of meta-analysis predictors
#'
#' @details
#' Test calculations for the predictors described in Table 52 of Environment Canada (2011) Scientific Assessment to Inform the Identification of Critical Habitat for Woodland Caribou (Rangifer tarandus caribou), Boreal Population, in Canada:2011 Update. Ottawa, Ontario.
#' So far, the predictors include:
#' \itemize{
#'   \item fire: % fire
#'   \item anthro: % non-overlapping anthropogenic disturbance.
#'   \item totalDist: Percent total non-overlapping fire and anthropogenic disturbance.
#' }
#'
#' Predictor values are calculated for a simple test disturbance data set that includes:
#'    fire: a single circular fire of radius = width/3
#'    anthroPoly: a single circular anthropogenic disturbance polygon of radius=width/3
#'    linear: a single circular road of radius=width
#'
#' Predictor values are returned for 3 buffer widths:
#'   answersNoBuffer - buffer width = 0
#'   answersBuffer - buffer width = width
#'   answersBufferHalf - buffer width = width/2
#'
#' @param width Numeric. Radius of circular linear
#' @param makeRaster vector of character strings. Optional. The set of maps to convert to rasters. Options include: ranges, anthroPoly, fires
#' @param nrows Numeric. Number of rows in output rasters. Ignored if makeRaster=c().
#' @param ncols Numeric. Number of columns in output rasters. Ignored if makeRaster=c().
#' @param forceRaster Boolean. Default F. Use for testing raster forceRaster option in landscapeStateBCaribou() method.
#' @param useLSClass Boolean. Default T. If T cLS is a landscapeStateBCaribou object. Otherwise cLS is a list of maps.
#' @return A list. cLS is a list of disturbance maps or landscapeStateBCaribou object. answersBuffer, answersNoBuffer, and answersBufferHalf contain the predictor values.
#' @examples
#'
#' require(raster)
#' cWidth=300
#' testData =makeTestPolygons(cWidth)
#' raster::plot(testData$cLS$ranges)
#' raster::plot(testData$cLS$fires,add=T,col="red")
#' raster::plot(testData$cLS$anthroPoly,add=T,col="blue")
#' raster::plot(testData$cLS$linear,add=T,col="black")
#' #meta-analysis predictors for buffer width = 300
#' testData$answersBuffer
#'
#' @export
makeTestPolygons<-function(width=300,makeRaster=c(),nrows=1000,ncols=1000,forceRaster=F,useLSClass=T){
  #Make a geometrically simple disturbance data set for testing.
  #width=300;makeRaster=c();nrows=1000;ncols=1000;forceRaster=F;useLSClass=T

  fireRadius =width/3
  anthroRadius = fireRadius
  roadRadius = width
  rangeRadius = 3*width

  basePointMap = sp::SpatialPoints(data.frame(x=0,y=0),proj4string=sp::CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  firePointMap = sp::SpatialPoints(data.frame(x=fireRadius+anthroRadius,y=0),proj4string=raster::crs(basePointMap))

  cRanges = raster::buffer(basePointMap,rangeRadius)
  if(length(makeRaster)>0){
    baseRaster = raster::raster(raster::extent(cRanges),nrows=1000,ncols=1000,crs=raster::crs(cRanges))
  }
  if(is.element("ranges", makeRaster)){
    cRanges = raster::rasterize(cRanges,baseRaster)
  }

  cAnthroPoly = raster::buffer(basePointMap,anthroRadius)
  if(is.element("anthroPoly", makeRaster)){
    cAnthroPoly =raster::rasterize(cAnthroPoly,baseRaster)
  }

  cFires = raster::buffer(firePointMap,fireRadius)
  if(is.element("fires", makeRaster)){
    cFires =raster::rasterize(cFires,baseRaster)
  }
  #plot(cLS@maps$ranges)
  #plot(cLS@maps$fires[[1]],add=T,col="red")
  #plot(cLS@maps$anthroPoly[[1]],add=T,col="blue")
  #anthroArea = cellStats(cLS@maps$bufferAnthro,stat="sum")
  #rangeArea = cellStats(cLS@maps$ranges,stat="sum")
  #fireArea =  cellStats(cLS@maps$fires,stat="sum")

  cLinear = raster::buffer(basePointMap, roadRadius)
  cLinear <- as(cLinear, 'SpatialLines')
  cForceRaster=forceRaster

  if(useLSClass){
    cLS = landscapeStateBCaribou(ranges=cRanges,fires=cFires,anthroPoly=cAnthroPoly,linear=cLinear,forceRaster=cForceRaster)
  }else{
    cLS = list(ranges=cRanges,fires=cFires,anthroPoly=cAnthroPoly,linear=cLinear)
  }
  rangeArea = pi*rangeRadius^2
  fireArea = pi*fireRadius^2

  answersNoBuffer = data.frame(fire=100*fireArea/rangeArea,anthro=100*pi*anthroRadius^2/rangeArea)
  answersNoBuffer$totalDist = answersNoBuffer$fire+answersNoBuffer$anthro

  answersBuffer = answersNoBuffer
  answersBuffer$anthro = 100*pi*(roadRadius+width)^2/rangeArea - answersBuffer$fire
  answersBuffer$totalDist = answersBuffer$anthro+answersBuffer$fire

  answersBufferHalf = answersBuffer
  answersBufferHalf$anthro = 100*pi*(roadRadius+width*0.5)^2/rangeArea - answersBuffer$fire
  answersBufferHalf$totalDist = answersBufferHalf$anthro+answersBufferHalf$fire

  return(list(cLS=cLS,answersBuffer=answersBuffer,answersNoBuffer=answersNoBuffer,answersBufferHalf=answersBufferHalf))
}

lsUnion<-function(mapA,mapB,mapAWarning="mapA",mapBWarning="mapB"){
  #mapA=rMaps[[cName]];mapB=lMaps[[cName]]
  ret=reconcileMapTypes(mapA,mapB,mapAWarning,mapBWarning)
  mapA=ret$mapA;mapB=ret$mapB

  if(checkNotEmpty(mapA)|checkNotEmpty(mapB)){
    if(ret$mapType=="raster"){
      mapA[is.na(mapA)]=0;mapB[is.na(mapB)]=0
      mapC = mapA|mapB
    }else{
      mapC=rgeos::gUnion(mapA,mapB)
    }
  }else{
    mapC=mapA
  }
  return(mapC)
}

lsDifference<-function(mapA,mapB,mapAWarning="mapA",mapBWarning="mapB"){
  #mapA=rMaps[[cName]];mapB=lMaps[[cName]]
  ret=reconcileMapTypes(mapA,mapB,mapAWarning,mapBWarning)
  mapA=ret$mapA;mapB=ret$mapB

  if(checkNotEmpty(mapB)&&checkNotEmpty(mapA)){
    if(ret$mapType=="raster"){
      #  if(is.element(1,raster::unique(x@maps$fires[[cName]]))){
      mapC=mapA
      mapC[mapB>0]=0
    }else{
      mapC=rgeos::gDifference(mapA,mapB)
    }
  }else{
    mapC=mapA
  }
  return(mapC)
}

checkNotEmpty<-function(map){
  #map=x[[kk]]
  if(is.element(class(map),c("RasterLayer","RasterBrick"))){
    return(length(setdiff(raster::unique(map),0))>0)
  }
  if(is.element(class(map),c("SpatialPolygons","SpatialPolygonsDataFrame"))){
    return(length(map@polygons)>0)
  }
  if(is.element(class(map),c("SpatialLines","SpatialLinesDataFrame"))){
    return(length(map@lines)>0)
  }
  stop("map class not recognized. Expecting Raster*, SpatialPolygons*, or SpatialLines* object.")
}

fillMissingMap<-function(cMaps=list(),cName,ranges){
  #ranges=x@maps$ranges
  if(!is.element(cName,names(cMaps))){
    if(is.element('static',names(cMaps))){
      cMaps[[cName]]=cMaps[['static']]
    }else{
      baseMap=ranges
      if(class(baseMap)=="RasterLayer"){
        cMaps[[cName]]=baseMap
        cMaps[[cName]][]=0
      }else{
        cMaps[[cName]] = baseMap[F,] #empty subset
      }
    }
  }
  return(cMaps)
}


.getFakePolyData<-function(inMap){
  ret=list(mapType="vector")
  #borrowing from https://gis.stackexchange.com/questions/140504/extracting-intersection-areas-in-r
  ret$mapB <- inMap
  # Remove attribute data
  ret$mapB <- as(ret$mapB, 'SpatialPolygons')
  # remove some of the polygons
  ret$mapB = ret$mapB[c(8,7,1,11,6,10,12,9),]

  # range polygons
  ret$mapA = raster::union(as(raster::extent(6, 6.4, 49.75, 50), 'SpatialPolygons'),
                   as(raster::extent(5.8, 6.2, 49.5, 49.7), 'SpatialPolygons'))
  #ret$mapA <- SpatialPolygonsDataFrame(ret$mapA, data.frame(field=c('x','y')), match.ID=F)
  raster::projection(ret$mapA) <- raster::projection(ret$mapB)

  return(ret)
}

reconcileMapTypes<-function(mapA,mapB,mapAWarning="mapA",mapBWarning="mapB"){
  #mapA =rMaps[[cName]];mapB=x@maps$fires
  #if mapA or mapB is a raster, force both maps to raster.
  if(is.element(class(mapA),c("RasterLayer","RasterBrick"))){
    cMapType = "raster"
  }else{
    cMapType= "vector"
  }

  if(is.element(class(mapB),c("RasterLayer","RasterBrick"))){
    fMapType = "raster"
  }else{
    fMapType= "vector"
  }

  #if there is a map type mismatch, force to raster
  if ((cMapType=="vector")&&(fMapType=="raster")){
    mapA=raster::rasterize(mapA,mapB)
    cMapType="raster"
    if(!is.null(mapAWarning)){warning("Map type mismatch. Converting ",mapAWarning," to raster.")}
  }
  if ((cMapType=="raster")&&(fMapType=="vector")){
    if(checkNotEmpty(mapB)){
      mapB=raster::rasterize(mapB,mapA)
    }else{
      mapB =mapA;mapB[]=NA
    }
    fMapType="raster"
    if(!is.null(mapBWarning)){warning("Map type mismatch. Converting ",mapBWarning," to raster.")}
  }

  return(list(mapA=mapA,mapB=mapB,mapType=cMapType))
}

#check range map validity
checkRanges<-function(maps){
  if(!is.element("ranges",names(maps))){
    stop("Expecting ranges among the maps.")
  }
  if(!is.element(class(maps$ranges),c('RasterLayer','SpatialPolygons','SpatialPolygonsDataFrame'))){
    stop("Range map must be a RasterLayer or SpatialPolygons* object.")
  }
}

#calculate percentage of disturbance within each range.
rangeCrosstab<-function(rangeMap,distMap,distName){
  #rangeMap=cRanges;distMap=landscapeState@maps$bufferAnthro[["t1"]];distName="anthro"
  #rangeMap=ret$mapA;distMap=ret$mapB[[1]];distName="fire"
  outNames =c("range",distName)
  if(is.element(class(rangeMap),c("RasterBrick","RasterLayer"))){

    distMap[is.na(distMap)]=0
    class(rangeMap)
    ct = as.data.frame(raster::crosstab(rangeMap,distMap)) #useNA argument appears to be not working
    names(ct)=c("Var1","Var2","Freq")
    ct=subset(ct,!is.na(Var1))

    ct$Var2=as.numeric(as.character(ct$Var2))

    rangeSums = plyr::ddply(ct,"Var1",plyr::summarize,rangeSum =sum(Freq))
    ct=merge(ct,rangeSums,all.y=T)
    ct$Freq[is.na(ct$Freq)]=0
    ct$pct = 100*ct$Freq/ct$rangeSum
    ct=subset(ct,Var2==1)
    ct = subset(ct,select=c(Var1,pct))
    names(ct)=outNames
  }else{
    ct=crosstabPolgyons(rangeMap,distMap,outNames)
  }
  return(ct)
}
#' @export
crosstabPolgyons<-function(mapA,mapB,outNames=c("range","area")){
  # Get percentage of each mapA polygon covered by mapB polygons.
  #ret=.getFakePolyData(shapefile(system.file("external/lux.shp", package="raster")));mapA=ret$mapA;mapB=ret$mapB[1,]
  #mapA=rangeMap;mapB=distMap

  if(class(mapA)!="SpatialPolygonsDataFrame"){
    dat = data.frame(rrrID=sapply(slot(mapA, "polygons"), function(x) slot(x, "ID")))
    rownames(dat)=dat$rrrID
    mapA = sp::SpatialPolygonsDataFrame(mapA, dat,match.ID=T)
  }else{
    mapA$rrrID=rownames(mapA@data)
  }
  if(class(mapB)=="SpatialPolygonsDataFrame"){
    mapB$rrrID=NULL;mapB$areaAll=NULL
  }

  if(length(mapB@polygons)==0){iMap=NULL}else{iMap=raster::intersect(mapB, mapA)}
  if(is.null(iMap)){
    cOut = data.frame(range=rownames(mapA@data),pct=0)
    names(cOut)=outNames
    return(cOut)
  }
  #plot(mapB, axes=T); plot(mapA, add=T); plot(iMap, add=T, col='red')

  # Extract areas from polygon objects then attach as attribute
  iMap$area = raster::area(iMap)
  mapA$areaAll = raster::area(mapA)
  iMap@data=merge(iMap@data,subset(mapA@data,select=c(rrrID,areaAll)),all.y=T)
  iMap$area[is.na(iMap$area)]=0
  cSums =plyr::ddply(iMap@data,plyr::.(rrrID),plyr::summarize,coveredArea=sum(area),totalArea=max(areaAll))

  cSums$pct = 100*cSums$coveredArea/cSums$totalArea
  cSums=subset(cSums,select=c(rrrID,pct))
  names(cSums)=outNames

  return(cSums)
}

#' @export
.adjustBrickNames<-function(x,ctag="X",cType="this",returnNamesOnly=F){
  #x=myTransitionGroup;ctag=tag
  cNames = names(x);cNames=gsub(ctag,"",cNames,fixed=T)

  outTag = "t"

  nNames =suppressWarnings(as.integer(cNames))

  if (is.element(NA,nNames)){
    if(length(nNames)==1){

      #if given a single unnamed map, assume it is static over time.
      warning(paste0("Name cannot be interpreted as an integer so assuming ",cType," doesn't change over time."))
      names(x) = "static"
    }else{

      #assume this is an ordered time-series and rename
      warning(paste0("Names cannot be interpreted as integers, so assuming ",cType," is an ordered time-series starting at time 1."))
      names(x) = paste0(outTag,seq(1,length(cNames)))
    }
    return(x)
  }

  outOrder = data.frame(outName = nNames,outID=seq(1:length(nNames)))
  outOrder = outOrder[order(outOrder$outName),]
  names(x)=paste0(outTag,nNames)

  if(returnNamesOnly){
    return(names(x))
  }

  if(identical(outOrder$outName,nNames)){
    return(x)
  }

  if(class(x)=="list"){
    return(x[outOrder$outName])

  }
  if(class(x)=="RasterStack"){
    x = raster::brick(x)
  }
  if(class(x)=="RasterBrick"){
    return(raster::subset(x,paste0(outTag,outOrder$outName)))
  }

  stop("Class of x not recognized. Expecting a list or RasterBrick")
}

