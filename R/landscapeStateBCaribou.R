#' @include AAAClassDefinitions.R
#' @include helperFns.R
NULL

# @name LandscapeStateBCaribou
# @rdname LandscapeStateBCaribou-class
setMethod(f = 'initialize', signature = "LandscapeStateBCaribou", definition = function(.Object, maps,params=list(),tag="X",forceRaster=F) {
  #maps = list(ranges=myStratum[[1]],fires=distMaps[["Fire"]],anthroPoly=distMaps[["Harvest"]]);tag="tg_.it1.ts"
  #maps = list(ranges=ranges,fires=fires,anthroPoly=anthroPoly,linear=linear);tag="tg_.it1.ts";forceRaster=F
  #maps=list(ranges=cLS@maps$ranges,fires=cLS@maps$fires,anthroPoly=cLS@maps$anthroPoly,linear=cLS@maps$linear);tag="static";forceRaster=T
  recognizedMaps = c("ranges","fires","anthroPoly","linear","newLandings")
  oddMaps = setdiff(names(maps),recognizedMaps)
  if(length(oddMaps)>0){
    stop("Unrecognized map types: ",paste(oddMaps,collapse=","),". Recognized map types include ",paste(recognizedMaps,collapse=","),".")
  }

  if(length(tag)==1){tag=list(fires=tag,anthroPoly=tag,linear=tag,newLandings=tag)}

  ret=checkRanges(maps)

  defNames =setdiff(recognizedMaps,"ranges")

  allTimes = c()
  for(cN in defNames){
    #cN="fires"

    if(!is.element(cN,names(maps))||is.null(maps[[cN]])){
      if(cN=="linear"){
        maps[[cN]] = SpatialLines(list(Lines(Line(cbind(x=c(0,1),y=c(0,1))), ID="a")),proj4string=crs(maps$ranges))
        maps[[cN]]=maps[[cN]][F,]
      }else{
        if(class(maps$ranges)=="RasterLayer"){
          maps[[cN]]=raster::brick(maps$ranges)
          maps[[cN]][]=0
          names(maps[[cN]])='fred'
        }else{
          maps[[cN]] = maps$ranges[F,] #empty subset
        }
      }
    }
    #allTimes=union(allTimes, gsub(tag[[cN]],"",names(maps[[cN]])))
  }
  #allTimes=setdiff(suppressWarnings(as.numeric(allTimes)),c(NA))

  ret = new("LandscapeState",maps=maps,tag=tag,forceRaster=forceRaster)
  .Object@maps=ret@maps
  .Object@params=ret@params
  .Object@mapType = ret@mapType
  return(.Object)
})


#' Landscape state for boreal caribou.
#' Returns a LandscapeStateBCaribou object containing correctly named and alligned maps for use by bCaribouMetaAnalysisPredictors().
#'
#' @details
#' All input maps must be alligned.
#' For raster fires and anthroPoly layers: NA and 0 areas will both be interpreted as no disturbance. Areas to be omitted from analysis (e.g. reservoirs) should be noted as NA in the range map. All non 0/NA values will be interpreted as disturbed (TRUE).
#' For ranges, fires, anthroPoly, and linear: if a single map is given, and names() cannot be interpreted as integers, will assume the feature is static over time.
#'
#' @param ranges Raster* or SpatialPolygon*. Range IDs. NA values will be excluded from analysis.
#' @param fires Raster*, SpatialPolygon*, or NULL. One or more maps indicating fire disturbance (T,F, or NA). If names cannot be interpreted as integers assume an ordered time-series. For rasters, NA and 0 areas will both be interpreted as no disturbance. Areas to be omitted from analysis (e.g. reservoirs) should be omitted from the range map.
#' @param anthroPoly Raster*, SpatialPolygon*, or NULL. One or more maps indicating polygonal anthropogenic disturbance (T,F, or NA). If names cannot be interpreted as integers assume an ordered time-series. See "tag" parameter for more details. For rasters, NA and 0 areas will both be interpreted as no disturbance. Areas to be omitted from analysis (e.g. reservoirs) should be omitted from the range map.
#' @param linear SpatialLines*, Raster* or NULL. One or more linear disturbance maps. If names cannot be interpreted as integers assume an ordered time-series. See "tag" parameter for more details. Areas to be omitted from analysis (e.g. reservoirs) should be omitted from the range map.
#' @param newLandings SpatialPoint*, SpatialPolygon*, Raster* or NULL. One or more maps of new features to be connected to the road network. If names cannot be interpreted as integers assume an ordered time-series. See "tag" parameter for more details. 
#' @param tag Character string or list of these. The part of names(fires) and names(anthroPoly) that is not the year. For example, if names(fires)=c("y3","y4") and names(anthroPoly)=c("t3","t4"), set tag=list(fires="y",anthroPoly="t").
#' @param forceRaster Boolean. Default F. If TRUE input maps will be converted to rasters alligned with ranges input map.
#' @return A validated LandscapeStateBCaribou object, containing alligned maps that are correctly named.
#' @examples
#' # TODO: examples
#' @export
setGeneric('landscapeStateBCaribou',function(ranges,fires=NULL,anthroPoly=NULL,linear=NULL,newLandings=NULL,tag="X",forceRaster=F) standardGeneric('landscapeStateBCaribou'))

#' @rdname landscapeStateBCaribou
setMethod('landscapeStateBCaribou', signature(ranges="RasterLayer"), function(ranges,fires,anthroPoly,linear,newLandings,tag,forceRaster) {
  #fires=NULL
  ret = new("LandscapeStateBCaribou",maps=list(ranges=ranges,fires=fires,anthroPoly=anthroPoly,linear=linear,newLandings=newLandings),tag=tag,forceRaster=forceRaster)
  return(ret)
})

#' @rdname landscapeStateBCaribou
setMethod('landscapeStateBCaribou', signature(ranges="SpatialPolygons"), function(ranges,fires,anthroPoly,linear,newLandings,tag,forceRaster) {
  #fires=NULL
  ret = new("LandscapeStateBCaribou",maps=list(ranges=ranges,fires=fires,anthroPoly=anthroPoly,linear=linear,newLandings=newLandings),tag=tag,forceRaster=forceRaster)
  return(ret)
})
