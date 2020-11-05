#' @include AAAClassDefinitions.R
#' @include helperFns.R
NULL

#' Apply buffer.
#'
#' @details
#' For a Raster* object use raster focal.
#'
#' For a LandscapeStateBCaribou* object, this function calculates non-overlapping buffered anthropogenic disturbance maps.
#' For more details, see appendix 7.5 of Environment Canada (2011) Scientific Assessment to Inform the Identification of Critical Habitat for Woodland Caribou (Rangifer tarandus caribou), Boreal Population, in Canada:2011 Update. Ottawa, Ontario.
#'
#' @param x Raster* or LandscapeStateBCaribou object.
#' @param width Numeric. Width of buffer to apply, in the units of x. Passed to raster::buffer.
#' @param subset Logical expression. Optional. Passed as an argument to the subset function in order to select a subset of raster values or polygons. Default for rasters is expression(x!=0) - 0s and NAs are ignored. Default for polygons is NULL (all polygons selected).
#' @examples
#' # TODO: examples
#' @export
setGeneric('lsBuffer',function(x,width,subset=NULL) standardGeneric('lsBuffer'))

#' @return A binary RasterLayer object showing buffered area.
#' @rdname lsBuffer
#' @export
setMethod('lsBuffer', signature(x="RasterLayer"), function(x,width,subset) {
  #Note there is a bug in raster::buffer - will be fixed in next version.
  #https://gis.stackexchange.com/questions/264133/raster-buffer-error-with-package-updates
  #x=testData$cLS@maps$anthroPoly[[1]];subset=NULL;width=300#x=x[[kk]]


  if(is.null(subset)){subset=expression(x!=0)}

  xtemp =x #distinguish between NA values and 0 values. NA values are omitted from calculation.
  x[is.na(x)]=0

  x[!eval(subset)]=0 #TO DO test this.

  if(length(raster::unique(x))>0){
    if(width>0){

      x=movingWindowAvg(x,radius=width,nms="name")#TO DO: check this works.
      #hdim = ceiling(width/res(x)[1])
      #wDim = 2*hdim+1
      #locSeq = seq(-hdim,hdim)
      #y = matrix(locSeq, wDim, wDim)
      #xx=t(y)
      #weights = (xx^2+y^2)^0.5
      #weights=weights <= width/res(x)[1]

      #vx = velox::velox(x)
      #vx$sumFocal(weights=weights, bands=c(1))

      #x = vx$as.RasterLayer()
      x[x!=0]=1
    }else{
      x[x>1]=1
    }
  }
  x[is.na(x)]=0
  x[is.na(xtemp)]=NA
  return(x)
})

#' @return A binary RasterBrick object.
#' @rdname lsBuffer
#' @export
setMethod('lsBuffer', signature(x="RasterStack"), function(x,width,subset) {
  #x=cLS@maps$anthroPoly
  return(lsBuffer(raster::brick(x),width,subset))
})

#' @return A binary RasterBrick object.
#' @rdname lsBuffer
#' @export
setMethod('lsBuffer', signature(x="RasterBrick"), function(x,width,subset) {
  #x=cLS@maps$anthroPoly;width=100;subset=NULL
  #checkNotEmpty<-function(x){T}
  #can't find a way to vectorize this
  xNames = names(x)
  for(kk in 1:raster::nlayers(x)){
    #kk=1
    if(checkNotEmpty(x[[kk]])){
      x[[kk]]=lsBuffer(x[[kk]],width=width,subset=subset)
    }
  }
  names(x)=xNames

  #out=raster::stackApply(x,indices=seq(1:nlayers(x)),fun=raster::buffer,width=width)
  return(x)
})

#' @return A SpatialPolygons object showing buffered area.
#' @rdname lsBuffer
#' @export
setMethod('lsBuffer', signature(x="SpatialPolygons"), function(x,width,subset) {
  #x=cLS@maps$anthroPoly[[1]];width=bufferWidth;subset=expression(Class=='Cutblock')

  if(!is.null(subset)){
    x=subset(x,eval(subset))
  }

  if(length(x@polygons)>0){
    x = raster::buffer(x,width=width)
  }
  return(x)
})

#' @return A SpatialPolygons object showing buffered area.
#' @rdname lsBuffer
#' @export
setMethod('lsBuffer', signature(x="SpatialLines"), function(x,width,subset) {
  #x=cLS@maps$linear[[1]];width=0;subset=expression(Class=='Road')


  if(!is.null(subset)){
    x=subset(x,eval(subset))
  }

  if(length(x@lines)>0){
    #Can run into trouble with large road networks.
    #"TopologyException: No forward edges found in buffer subgraph"
    #https://gis.stackexchange.com/questions/31880/memory-issue-when-trying-to-buffer-union-large-dataset-using-postgis

    if(width==0){
      x = rgeos::gConvexHull(x)
      x=x[F,]
    }else{
      temp = raster::buffer(x,width=width,dissolve=F)
      x =rgeos::gUnaryUnion(temp)
    }
  }
  return(x)
})

#' @return A list of buffered objects.
#' @rdname lsBuffer
#' @export
setMethod('lsBuffer', signature(x="list"), function(x,width,subset) {
  for(i in 1:length(x)){
    x[[i]]=lsBuffer(x[[i]],width=width,subset=subset)
  }
  return(x)
})

#' @return A LandscapeStateBCaribou object with filled anthroBuffer slot. Overlapping fire disturbed areas are removed from the buffered area.
#' @rdname lsBuffer
#' @export
setMethod('lsBuffer', signature(x="LandscapeStateBCaribou"), function(x,width,subset) {
  #x=cLS;subset=NULL;width=500

  #Buffer layers
  print("Buffering anthropogenic disturbance polygons.")
  rMaps = lsBuffer(x@maps$anthroPoly,width=width,subset=subset)

  #rMaps=list()
  if(is.element("linear",names(x@maps))){
    print("Buffering linear features.")
    lMaps = lsBuffer(x@maps$linear,width=width,subset=subset)
  }else{
    lMaps =list()
  }

  allNames=union(names(rMaps),names(lMaps))
  if(!identical('static',allNames)){
    allNames=setdiff(allNames,'static')
  }
  #remove areas that are also disturbed by fire, and merge buffered areas

  for(cName in allNames){
    #cName="t4"

    #add empty maps if necessary
    rMaps=fillMissingMap(rMaps,cName,x@maps$ranges)
    lMaps=fillMissingMap(lMaps,cName,x@maps$ranges)
    cFires=fillMissingMap(x@maps$fires,cName,x@maps$ranges)

    rMaps[[cName]]=lsUnion(rMaps[[cName]],lMaps[[cName]],
                           paste0("buffered anthroPoly map",cName),
                           mapBWarning=NULL) #expect linear disturbance conversion here so don't warn of mismatch.
    rMaps[[cName]]=lsDifference(rMaps[[cName]],cFires[[cName]],
                                  paste0("buffered anthroPoly map",cName),
                                  paste0("fire map",cName))
    #Don't accumulate map in memory when they are no longer needed.
    if(class(cFires)=="list"){
      cFires[[cName]]=NULL
    }else{
      cFires=raster::dropLayer(cFires,cName)
    }
    #lMaps[[cName]]=NULL
  }

  x@maps$bufferAnthro = rMaps
  x@bufferWidths$bufferAnthro = width
  x@params$subset=subset
  return(x)
})
