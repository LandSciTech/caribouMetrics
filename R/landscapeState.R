#' @include AAAClassDefinitions.R
#' @include helperFns.R
NULL

# @name LandscapeState
# @rdname LandscapeState-class
setMethod(f = 'initialize', signature = "LandscapeState", definition = function(.Object, maps,params=list(),tag="X",forceRaster=F) {
  #maps = list(ranges=ranges,fires=fires,anthroPoly=anthroPoly,linear=linear);tag="tg_.it1.ts";forceRaster=F

  typeCounts = list(raster=0,polygon=0)
  acceptedInputMapFormats = c('RasterLayer','RasterStack','RasterBrick','SpatialPolygons','SpatialPolygonsDataFrame','SpatialLines','SpatialLinesDataFrame')

  ret=checkRanges(maps)
  if(forceRaster&&(class(maps$ranges)!='RasterLayer')){
    stop("Range map must be a RasterLayer in order to convert polygons to rasters.")
  }

  if(is.element(class(maps$ranges),c("RasterBrick","RasterLayer"))){
    if(raster::nlayers(maps$ranges)>1){
      stop("Expecting a single range map.")
    }
    maps$ranges=maps$ranges[[1]]
  }

  stateMaps = setdiff(names(maps),c("ranges"))

  if(length(stateMaps)==0){
    stop("Expecting maps to define the landscape state, not only a range map.")
  }

  if(length(tag)==1){
    oTag = tag;tag=list()
    for(i in 1:length(stateMaps)){
      tag[[stateMaps[[i]]]]=oTag
    }
  }

  if(class(maps$ranges)=='RasterLayer'){
    rasterRanges = maps$ranges
    typeCounts$raster=typeCounts$raster+1
  }else{
    rasterRanges = NA
    typeCounts$polygon=typeCounts$polygon+1
  }

  #check that input maps are valid and alligned
  for(i in 1:length(stateMaps)){
    #i =1
    cM = stateMaps[[i]]

    #Convert Raster* to RasterBricks - will complain if these are not alligned Raster* objects.
    if(is.element(class(maps[[cM]]),c('RasterLayer','RasterStack'))){
      mNames = names(maps[[cM]])
      maps[[cM]]=raster::brick(maps[[cM]])
      names(maps[[cM]])=mNames
    }

    if(!is.element(class(maps[[cM]]),c('RasterBrick','list'))){
      maps[[cM]]=list(static=maps[[cM]])
      tag[[cM]]="static"
    }

    if(!is.element(class(maps[[cM]][[1]]),acceptedInputMapFormats)){
        stop("Format of ", cM," map not recognized. Accepted formats are: ",paste(acceptedInputMapFormats,collapse=","))
    }

    if((cM=="linear")&&!is.element(class(maps[[cM]][[1]]),c("SpatialLines","SpatialLinesDataFrame","RasterLayer"))){
      #maps=list(linear=outRoads);cM="linear"
      if(class(maps[[cM]][[1]])=="RasterLayer"){
        warning("Converting linear rasters to line segments. This will work for output from roads::projectRoads(), and raster::rasterize(SpatialLines), but not for rasters in general.")
        outSet =names(maps[[cM]])
        outMaps=list()
        for(ll in 1:length(outSet)){
          #ll=2
          cName =outSet[ll]

          cMap = maps[[cM]][[cName]]

          cLines = convertRasterToLineSegments(cMap)#roads::convertRasterToLineSegments(cMap)#remove dependency on roads package

          outMaps[[cName]]=cLines
        }
        maps[[cM]]=outMaps
      }else{
        stop("Expecting linear maps to be SpatialLines, SpatialLinesDataFrame, or Raster* objects. Not ",class(maps[[cM]]),".")
      }
    }

    if(forceRaster){
      if(is.element(class(maps[[cM]][[1]]),c("SpatialPolygons","SpatialPolygonsDataFrame"))){
        for(kk in 1:length(maps[[cM]])){
          maps[[cM]][[kk]]=raster::rasterize(maps[[cM]][[kk]],rasterRanges)
        }
        maps[[cM]]=raster::brick(maps[[cM]])
      }
    }

    if(class(maps[[cM]])=='RasterBrick'){
      typeCounts$raster=typeCounts$raster+1

      #Now check alignment - will complain if not alligned
      if(class(rasterRanges)!="RasterLayer"){
        rasterRanges = raster::rasterize(maps$ranges,maps[[cM]][[1]])
      }
      checkAllign = raster::compareRaster(maps$ranges,maps[[cM]][[1]])
      if(!checkAllign){
        stop("Problem with ",cM,". All input rasters must have the same same extent, number of rows and columns,
projection, resolution, and origin.")
      }

      #If names are character integers, make sure they are ordered. Otherwise, assume an ordered time-series and rename accordingly.
      maps[[cM]]=.adjustBrickNames(maps[[cM]],tag[[cM]],cType=cM)
    }

    if(class(maps[[cM]])=='list'){
      if(is.element(class(maps[[cM]][[1]]),c('SpatialPolygons','SpatialPolygonsDataFrame'))){
        typeCounts$polygon=typeCounts$polygon+1
      }

      for(kk in 1:length(maps[[cM]])){
        ret=tryCatch(rgeos::gIsValid(maps[[cM]][[kk]]), error = function(e) stop("Geometry error. ",cM,"[[",kk,"]] is not valid. ",e))
        #Check allignment. Complain if extent and crs do not match.
        #if(raster::extent(maps$ranges)!=raster::extent(maps[[cM]])){
        #  stop("Problem with ",cM,". All input maps must have the same extent.")
        #}
        if(!compareCRS(maps$ranges,maps[[cM]][[kk]])){
          stop("Problem with ",cM,"[[",kk,"]]. All input maps must have the same coordinate reference system.")
        }

      }
      maps[[cM]]=.adjustBrickNames(maps[[cM]],tag[[cM]],cType=cM) #sort the list
    }
  }
  #Now assume maps are ordered RasterBricks or lists of Spatial* objects, ranges is a RasterLayer or SpatialPolygons* object, and everying is alligned.
  .Object@maps=maps
  .Object@params=params
  if((typeCounts$raster>0)&(typeCounts$polygon>0)){
    .Object@mapType='mixed'
  }else{
    if(typeCounts$raster>0){
      .Object@mapType="raster"
    }else{
      .Object@mapType="polygon"
    }
  }

  return(.Object)
})

