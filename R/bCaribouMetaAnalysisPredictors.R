#' @include AAAClassDefinitions.R
#' @include helperFns.R
NULL

#' Calculate national meta-analysis predictors for boreal caribou
#'
#' @details
#' Calculate the predictors described in Table 52 of Environment Canada (2011) Scientific Assessment to Inform the Identification of Critical Habitat for Woodland Caribou (Rangifer tarandus caribou), Boreal Population, in Canada:2011 Update. Ottawa, Ontario.
#' So far, the variables calculated by this function include:
#' \itemize{
#'   \item fire: % fire
#'   \item anthro: % non-overlapping anthropogenic disturbance. Buffer in landscapeState@maps$bufferWidth.
#'   \item totalDist: Percent total non-overlapping fire and anthropogenic disturbance.
#' }
#'
#' Note that NA values are omitted from tabulated area.
#' Missing layers are omitted from the output, not interpreted as 0 disturbance.
#'
#' @param landscapeState LandscapeStateBCaribou object. Contains the necessary range and disturbance maps.
#' @param years Integer or vector of these. Default is all the years included in landscapeState object.
#' @param ranges Integer, character or vector of these. IDs of ranges to include in the analysis. Default is all the ranges in landscapeState object.
#' @param predictors Character string or vector of these. Output variables to calculate. See details for options.
#' @param width Numeric. Optional. Width of buffer to apply - see lsBuffer for details. If NULL (default) buffer information within landscapeState object will be used if available. If width is specified buffer information within the landscapeState object will be ignored.
#' @return A data frame of predictor variables for each year and range.
#' @examples
#' # TODO: examples
#' @export
setGeneric('bCaribouMetaAnalysisPredictors',function(landscapeState,years=NULL,ranges=NULL,predictors=c("fire","anthro","totalDist"),width=NULL) standardGeneric('bCaribouMetaAnalysisPredictors'))

#' @rdname bCaribouMetaAnalysisPredictors
setMethod('bCaribouMetaAnalysisPredictors', signature(landscapeState="LandscapeStateBCaribou"), function(landscapeState,years,ranges,predictors,width) {
  #landscapeState = cLS;years=NULL;ranges=NULL;predictors=c("anthro");width=300

  recognizedPredictors = c("fire","anthro","totalDist")
  bufferPredictors = c("anthro","totalDist")


  #if width is specified, redo buffer analysis
  if(!is.null(width)){
    landscapeState = lsBuffer(landscapeState,width=width)
  }


  if((length(intersect(bufferPredictors,predictors))>0)&&!is.element("bufferAnthro",names(landscapeState@maps))){
    stop("Buffering not yet done. Please specify width parameter or call lsBuffer() before calling bcMetaAnalysisPredictors().")
  }

  #disturbances are already non-overlapping. So simply a matter of summarizing area by range, and adding.

  #Get years
  allYrs = getAllYrs(union(names(landscapeState@maps$fires),names(landscapeState@maps$anthroPoly)))

  if(!is.null(years)){
    #years = seq(90:100)
    missingYrs = setdiff(years,allYrs)
    if(length(missingYrs)>0){
      warning("landscapeState does not contain disturbance data for selected years: ",paste(missingYrs,sep=","))
    }
    includeYrs = sort(intersect(years,allYrs))
    if(length(includeYrs)==0){
      stop("Nothing to calculate. landscapeState does not contain disturbance data for any of the selected years.")
    }
  }else{
    includeYrs = allYrs
  }

  #get revised range map
  cRanges = landscapeState@maps$ranges

  if(class(cRanges)=="RasterLayer"){
    allRanges = raster::unique(cRanges)
  }else{
    allRanges = sapply(slot(cRanges, "polygons"), function(x) slot(x, "ID"))
    if(!is.null(ranges)){ranges = as.character(ranges)}
  }

  if(!is.null(ranges)){
    #ranges = c(1)

    missingRanges = setdiff(ranges,allRanges)
    if(length(missingRanges)>0){
      warning("landscapeState@maps$ranges does not contain selected ranges: ",paste(missingRanges,sep=","))
    }
    includeRanges = intersect(ranges,allRanges)
    if(length(includeRanges)==0){
      stop("Nothing to calculate. landscapeState@maps$ranges does not contain any of the selected ranges.")
    }
    if(class(cRanges)=="RasterBrick"){
      df = data.frame(by=includeRanges,which=includeRanges)
      cRanges = raster::subs(cRanges,df)
    }else{
      cRanges = cRanges[includeRanges,]
    }
  }else{
    includeRanges=allRanges
  }

  #Check that selected predictors are valid
  missingPreds = setdiff(predictors,recognizedPredictors)
  if(length(missingPreds)>0){
    stop("Selected predictors not recognized: ",paste(missingPreds,sep=","))
  }

  #Now for each year and range get summary info
  outInfo = NULL
  for (y in includeYrs){
    #y=includeYrs[1]#'static'
    if(y!="static"){cMapTag = paste0("t",y)}else{cMapTag="static"}
    cSummary = data.frame(year=y)
    if(length(intersect(c("fire","totalDist"),predictors))>0){
      presentMap = intersect(c(cMapTag,"static"),names(landscapeState@maps$fires))
      if(length(presentMap)>0){
        cFire = landscapeState@maps$fires[[presentMap]]
        #plot(landscapeState@maps$fires)
        ret=reconcileMapTypes(cRanges, cFire)
        #ret=.getFakePolyData(shapefile(system.file("external/lux.shp", package="raster")));ret$mapB=ret$mapB[5,]
        #plot(ret$mapA)
        ct=rangeCrosstab(ret$mapA, ret$mapB,"fire")
        #plot(ret$mapB, axes=T); plot(ret$mapA, add=T)
      }else{
        ct = data.frame(range=includeRanges,fire=0)
      }
      cSummary=merge(cSummary,ct)
    }

    if(length(intersect(c("anthro","totalDist"),predictors))>0){
      #this is non-overlapping with fire. Fire areas are removed in the buffering step.
      presentMap = intersect(c(cMapTag,"static"),names(landscapeState@maps$bufferAnthro))
      if(length(presentMap)>0){
        ct=rangeCrosstab(cRanges, landscapeState@maps$bufferAnthro[[presentMap]],"anthro")
      }else{
        ct = data.frame(range=includeRanges,anthro=0)
      }
      cSummary=merge(cSummary,ct)

    }

    if(length(intersect(c("totalDist"),predictors))>0){
      present = intersect(c("fire","anthro"),names(cSummary))
      if(length(present)==2){
        cSummary$totalDist=cSummary$fire+cSummary$anthro
      }else{
        if(is.element("fire",present)){
          cSummary$totalDist=cSummary$fire
        }
        if(is.element("anthro",present)){
          cSummary$totalDist=cSummary$anthro
        }
      }
    }

    if(is.null(outInfo)){outInfo=cSummary}else{outInfo=rbind(outInfo,cSummary)}
  }

  return(outInfo)
})
