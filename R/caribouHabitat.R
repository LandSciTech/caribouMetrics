#' @include AAAClassDefinitions.R
NULL

setMethod(f = "initialize", signature = "CaribouHabitat",
          definition = function(.Object, landCover, esker, 
                                natDist, anthroDist, 
                                linFeat, projectPoly,  
                                processedData, habitatUse, attributes){
            .Object@landCover <- landCover 
            .Object@esker <- esker
            .Object@natDist <- natDist
            .Object@anthroDist <- anthroDist
            .Object@linFeat <- linFeat
            .Object@projectPoly <- projectPoly
            .Object@processedData <- processedData
            .Object@habitatUse <- habitatUse 
            .Object@attributes <- attributes
            return(.Object)
          })

#'caribouHabitat
#'
#'Calculate the probability of caribou habitat use in spring, summer, fall and
#'winter for caribou ranges in Northern Ontario, based on Hornseth and Rempel,
#'2016.
#'
#'Caribou habitat use is calculated based on the availability of resources and
#'the presence of disturbances on the landscape. The primary source of resource
#'information is the \code{landCover} but this is updated based on disturbance
#'information. All data sources can be provided either as filenames or as
#'spatial files. If filenames are provided then the \code{landCover} is assumed
#'to be the Provincial Landcover for Ontario and is converted to resource types
#'using \code{reclassPLC}. The result is a CaribouHabitat object which has
#'methods defined for plotting and extracting the results. To update an existing
#'CaribouHabitat object with new data see \link[caribouMetrics]{updateCaribou}.
#'
#'
#'@param landCover filename or RasterLayer. Provincial landcover class
#'@param esker filename, RasterLayer or sf object. Eskers. If it is a
#'  RasterLayer then it should be esker density in m^2/ha.
#'@param linFeat filename, RasterLayer, sf object or a list of these that will
#'  be combined. Linear features. If it is a RasterLayer then it should be
#'  linear feature density in m^2/ha. If a RasterLayer is provided as a list
#'  element then ptDensity will be used to assign a density of linear features
#'  in the pixel (default is 1).
#'@param projectPoly filename or sf object. Polygon defining the project area.
#'  If caribouRange is a data.frame this must have a column called Range with
#'  the name of the caribou range represented by the polygon which corresponds
#'  to the Range column in the caribouRange data.frame
#'@param caribouRange character or data.frame. If character the range where
#'  caribou were located. See \code{unique(coefTableHR$Range)} for options. If
#'  data.frame it must have two columns Range and coefRange. Range is the name
#'  of the geographical area and is used to link the table the the provided
#'  \code{projectPoly} polygons. coefRange is the name of the caribou range that
#'  the coefficients should be used from.
#'@param coefTable data.frame. table of coefficients to be used in the
#'   model. Must match the format and naming of the default \code{coefTableHR}
#'@param ... optional arguments:
#' \describe{
#'   \item{natDist:}{filename or RasterLayer. Presence or absence of natural 
#'   disturbance, primarily by fire.}
#'   \item{anthroDist:}{filename or RasterLayer. Anthropogenic disturbance including 
#'   harvest. This can have an effect on any type of land cover except water.}
#'   \item{eskerSave}{filename to save rasterized esker data.}
#'   \item{linFeatSave}{filename to save rasterized linear feature data.}
#'   \item{padProjPoly}{logical. Should the area around the \code{projectPoly} be
#'   used to avoid edge effects? If FALSE, the default, only data from inside the
#'   \code{projectPoly} is used. If TRUE then \code{projectPoly} is buffered and
#'   the other variables are clipped to the extent of the buffered area. Results
#'   are always clipped to the original \code{projectPoly}. It is ideal to set
#'   this to TRUE and provide a data set that is larger than the
#'   \code{projectPoly} to avoid edge effects.}
#'   \item{padFocal}{logical. This value is passed to the pad argument in
#'   \code{raster::focal}, if it is FALSE then cells near the edge will return
#'   NA, if it is TRUE a value will be returned for each cell that assumes cells
#'   outside the input data are 0 for all resource types. This is not a good
#'   assumption and should be used with caution.}
#'   \item{saveOutput}{character. The filename to save the rasterBrick of habitat
#'   use probabilities to. Note this will overwrite any existing files. The .grd
#'   format is recommended because it will preserve layer names when the file is
#'   reloaded.}
#'   \item{winArea}{number. This is the area of the moving window that is used to
#'   average proportions of each resource type at broader spatial scales. The
#'   Hornseth and Rempel (2016) models used specific window areas which are
#'   defined within this package and used as the default. You should only specify
#'   a window size if you have good reason.}
#'   \item{doScale}{logical. FALSE by default. Set to TRUE only if you have
#'   supplied coefficients that were trained on standardized data which will
#'   cause the input data to be scaled.}
#'   \item{ptDensity}{number. Only used if a list element in linFeat is a raster.
#'   See \code{\link{rasterizeLineDensity}}.}
#' }
#' 
#'@return A CaribouHabitat Object see \code{\link{CaribouHabitat-class}}
#'
#'@seealso \code{\link{CaribouHabitat-class}} for information on the object
#'  returned and \code{\link{updateCaribou}} for updating an existing
#'  CaribouHabitat object.
#'
#' @source Rempel, R.S., Carlson, M., Rodgers, A.R., Shuter, J.L., Farrell,
#'   C.E., Cairns, D., Stelfox, B., Hunt, L.M., Mackereth, R.W. and Jackson,
#'   J.M., 2021. Modeling Cumulative Effects of Climate and Development on
#'   Moose, Wolf, and Caribou Populations. The Journal of Wildlife Management.
#'
#'  Hornseth, M.L. and Rempel, R.S., 2016. Seasonal resource selection of
#'  woodland caribou (Rangifer tarandus caribou) across a gradient of
#'  anthropogenic disturbance. Canadian Journal of Zoology, 94(2), pp.79-93.
#'  \url{https://doi.org/10.1139/cjz-2015-0101}
#'
#'@export
setGeneric("caribouHabitat", 
           function(landCover, esker, linFeat, projectPoly, caribouRange, ...) 
             standardGeneric("caribouHabitat"))
#' @rdname caribouHabitat
setMethod(
  "caribouHabitat", 
  signature(landCover = "ANY"), 
  function(landCover, esker, linFeat, projectPoly, caribouRange, 
           coefTable = coefTableHR, ...) {

    dots <- list(...)
    
    # make sure caribouRange is a dataframe
    if(!inherits(caribouRange, "data.frame")){
      caribouRange <- data.frame(Range = caribouRange, 
                                 coefRange = caribouRange, 
                                 stringsAsFactors = FALSE)
    } else {
      if(any(names(caribouRange) != c("Range", "coefRange"))){
        stop("If caribouRange is a data.frame the column names",
             " must be Range and coefRange", call. = FALSE)
      }
    }
    
    # Get window area from table b/c some models used different sizes
    if(is.null(dots$winArea)){
      dots$winArea <- coefTable %>% filter(Range %in% caribouRange$coefRange) %>% 
        pull(WinArea) %>% 
        unique()
    }
    
    # If multiple winAreas need to apply processData separately
    if(length(dots$winArea) > 1){

      # reassign to NULL so it will be redetermined from each call
      dots$winArea <- NULL
      
      # provide template so resampled rasters will match 
      
      # NOTE: to adjust res of final model change this res and the res in
      # default argument in inputData
      dots$tmplt <- raster(landCover) %>% raster::`res<-`(c(400, 400))
      
      # polygons of ranges split into list with different winAreas
      projPolyLst <- projectPoly %>% 
        left_join(caribouRange, by = "Range") %>% 
        left_join(coefTable %>% group_by(Range) %>%
                    summarize(WinArea = first(WinArea)),
                  by = c(coefRange = "Range")) %>% 
        select(-coefRange) %>% 
        split(.$WinArea) %>% 
        purrr::map(~select(.x, -WinArea))
      
      # caribouRange values for each winArea
      carRangeLst <- caribouRange %>%
        left_join(coefTable %>% group_by(Range) %>%
                    summarize(WinArea = first(WinArea)),
                  by = c(coefRange = "Range")) %>% 
        split(.$WinArea) %>% 
        purrr::map(~select(.x, -WinArea))
   
      resultLst <- purrr::map2(projPolyLst, carRangeLst,
                        ~do.call(caribouHabitat, 
                                c(list(landCover = landCover, 
                                  esker = esker, 
                                  linFeat = linFeat, 
                                  projectPoly = .x, 
                                  caribouRange = .y, 
                                  coefTable = coefTable),
                                  dots)))
      
      # Re-combine CarHab objects into one
      x <- resultLst[[1]]
    
      # names of slots to iterate over
      slotLst <- slotNames(resultLst[[1]])
      slotLst <- slotLst[!slotLst %in% c("attributes", "projectPoly")]
      
      for(i in 1:length(slotLst)){
        slotNm <- slotLst[i]
        rastLst <- lapply(resultLst, slot, name = slotNm)
        slot(x, slotNm) <- doMosaic(rastLst)
      }
      
      # Only thing that will change about projectPoly is crs 
      x@projectPoly <- projectPoly %>% st_transform(st_crs(landCover))
      
    } else {
      
      inputDataArgs <- dots[c("natDist", "anthroDist", 
                              "winArea", "eskerSave", "linFeatSave", 
                              "padProjPoly", "padFocal", "ptDensity", 
                              "tmplt")]
      
      inputDataArgs <- inputDataArgs[which(lapply(inputDataArgs, length) > 0)]
      
      x <- do.call(inputData, c(lst(landCover, esker, linFeat, projectPoly,
                                    caribouRange, coefTable), 
                                inputDataArgs))
      
      x <- processData(x)
      
      updateArgs <- dots[c("coefTable", "doScale")]
      
      updateArgs <- updateArgs[which(lapply(updateArgs, length) > 0)]
      
      x <- do.call(updateCaribou, c(list(CarHab = x), updateArgs))
    }
    
    if(!is.null(dots$saveOutput)){
      
      byLayer <- grepl("\\.asc$|\\.sdat$|\\.rst$", dots$saveOutput)
      if(!byLayer && !grepl("\\.grd", dots$saveOutput)){
        warning("Saving output to ", dots$saveOutput, 
                ". Layernames will not be preserved.",
                " Use .grd format to preserve names")
      }
 
      raster::writeRaster(x@habitatUse, filename = dots$saveOutput, 
                          overwrite = TRUE, bylayer = byLayer, 
                          suffix = "names")
    }
    
    return(x)
  })

doMosaic <- function(rastLst){
  # do.call doesn't work with names
  names(rastLst) <- NULL
  
  rastLst$fun <- mean
  
  out <- do.call(raster::mosaic, rastLst)
  
  names(out) <- names(rastLst[[1]])
  
  return(out)
}