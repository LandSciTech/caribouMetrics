#' @include AAAClassDefinitions.R
NULL

#' Update an existing CaribouHabitat Object
#'
#' Update a CaribouHabitat object in order to avoid reprocessing parts of the
#' data that have not changed. New data is supplied as a named list and the
#' object is updated depending on the elements provided in the list.
#'
#' If `newData` contains only linFeat then only linear features will be
#' updated.
#'
#' If  `newData` contains landCover the landCover in the CaribouHabitat
#' object will be replaced and new projections made for the whole landscape. If
#' natDist or anthroDist are not provided then the data stored in the
#' CaribouHabitat object is reused.
#'
#' @param CarHab CaribouHabitat object
#' @param newData a named list of objects to be used to update
#'   CarHab. Potential names are: landCover, natDist, anthroDist, and
#'   linFeat.
#' @param resultsOnly logical. If FALSE the whole CaribouHabitat object is
#'   returned. If TRUE then only the habitatUse RasterStack is returned.
#' @param coefTable data.frame. Optional table of coefficients to be used in the
#'   model. Must match the format and naming of `coefTableHR`
#' @param doScale logical. FALSE by default. Set to TRUE only if you have
#'   supplied coefficients that were trained on standardized data which will
#'   cause the input data to be scaled.
#' @param ... other arguments passed to methods. Not currently used.
#'
#' @return If `resultsOnly` is FALSE an updated CaribouHabitat object. If
#'   `resultsOnly` is TRUE a RasterStack with a layer for each season.
#'   
#' @examples 
#' # create example rasters
#' lc <- raster::raster(xmn = 0, xmx = 25000, ymn = 0, ymx = 25000, 
#'                      resolution = 250, crs = 5070)
#' lc[] <- 0
#' nd <- lc
#' nd[1:30, 1:30] <- 1
#' ad <- lc
#' ad[30:50, 3:50] <- 1
#' lc[] <- 1
#' lc[70:100, 70:100] <- 2
#' 
#' # create sf objects
#' lf <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10000, 10000),
#'                                                             ncol = 2, byrow = TRUE))),
#'                               crs = 5070))
#' esk <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 10000, 10000, 0),
#'                                                             ncol = 2, byrow = TRUE))),
#'                               crs = 5070))
#' 
#' 
#' projPol <- sf::st_sf(sf::st_as_sfc(sf::st_bbox(ad)))
#' 
#' # calculate relative probability of use
#' res <- caribouHabitat(landCover = lc,
#'                linFeat = lf,
#'                esker = esk,
#'                natDist = nd,
#'                anthroDist = ad,
#'                projectPoly = projPol,
#'                caribouRange = "Nipigon",
#'                winArea = 1000 #leave as default NULL except for small examples
#' )
#' 
#' # new linear features
#' lf2 <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 25000, 25000),
#'                                                              ncol = 2, byrow = TRUE)),
#'                                     sf::st_linestring(matrix(c(0, 20000, 20000, 25000),
#'                                                              ncol = 2, byrow = TRUE)),
#'                                     sf::st_linestring(matrix(c(10000, 0, 20000, 20000),
#'                                                              ncol = 2, byrow = TRUE))),
#'                               crs = 5070))
#' 
#' 
#' res2 <- updateCaribou(res, newData = list(linFeat = lf2))
#' 
#' # visualize the impact of new roads
#' plot(res, season = "Winter")
#' plot(res2, season = "Winter")
#' 
#' @family habitat
#' @export
setGeneric("updateCaribou", function(CarHab, newData, ...) standardGeneric("updateCaribou"))

#' @rdname updateCaribou
# method to get final calculation from processed data in CaribouHabitat object
setMethod(
  "updateCaribou", 
  signature(CarHab = "CaribouHabitat", newData = "missing"), 
  function(CarHab, coefTable = coefTableHR, doScale = FALSE) {
    # process the data if not already done
    if(nrow(CarHab@processedData) < 2){
      
      CarHab <- processData(CarHab)
    }
    

    # which coefficients to use for which range
    rangeCoefLst <- CarHab@projectPoly %>% 
      left_join(CarHab@attributes$caribouRange, by = "Range")
    
    if(length(unique(rangeCoefLst$coefRange)) > 1){
      
      rangeCoefLst <- split(rangeCoefLst, rangeCoefLst$coefRange)
      
      applyCalcRSP <- function(dat, rangeCoef, doScale, coefTable){
        dat <- terra::mask(dat, rangeCoef)
        
        coefT <- coefTable %>% 
          filter(.data$Range %in% rangeCoef$coefRange)
        
        habitatUse <- calcRSP(dat, coefT, doScale = doScale)
      }
      
      habUseLst <- lapply(rangeCoefLst, applyCalcRSP, dat = CarHab@processedData,
                          doScale = doScale, coefTable = coefTable)
    
      # do.call doesn't work with names
      names(habUseLst) <- NULL
      
      habUseLst$fun <- function(...){sum(..., na.rm = TRUE)}
      
      CarHab@habitatUse <- do.call(terra::lapp, habUseLst)

      # 0 created by sum when all are NA but reintroduce NA from processed
      getNA <- !is.na(CarHab@processedData$CON)
      getNA[getNA == 0] <- NA
      CarHab@habitatUse <- CarHab@habitatUse * getNA
      
      names(CarHab@habitatUse) <- names(habUseLst[[1]])
    } else {
      coefT <- coefTable %>% 
        filter(.data$Range %in% rangeCoefLst$coefRange)
      
      CarHab@habitatUse <- calcRSP(CarHab@processedData, coefT, 
                                   doScale = doScale)
       
    }

    # Takes a long time not sure it is worth it
    # # This keeps cells that are only partially in the polygon 
    # projRas <- raster::rasterize(CarHab@projectPoly, CarHab@habitatUse[[1]],
    #                              getCover=TRUE)
    # projRas[projRas==0] <- NA
    # 
    CarHab@habitatUse <- terra::crop(CarHab@habitatUse, CarHab@projectPoly)
    # CarHab@habitatUse <- raster::crop(CarHab@habitatUse, CarHab@projectPoly,
    #                                   snap = "out")
    # 
    CarHab@processedData <- terra::mask(CarHab@processedData, CarHab@projectPoly)
    # CarHab@processedData <- raster::crop(CarHab@processedData, CarHab@projectPoly,
    #                                      snap = "out")
    
    return(CarHab)
  })

# method to update processed data when new data is supplied 
#' @rdname updateCaribou
setMethod(
  "updateCaribou", 
  signature(CarHab = "CaribouHabitat", newData = "list"), 
  function(CarHab, newData, resultsOnly = FALSE, 
           coefTable = coefTableHR, doScale = FALSE) {

    if(nrow(CarHab@processedData) < 2){
      stop("CarHab@processedData is empty. Run updateCaribou with no additional
           data to process the initial data before updating")
    }
    
    
    CarHab <- processData(CarHab, newData)
    
    CarHab <- updateCaribou(CarHab, coefTable = coefTable, doScale = doScale)
    
    if(resultsOnly){
      return(CarHab@habitatUse)
    }
    
    return(CarHab)
  
  })

