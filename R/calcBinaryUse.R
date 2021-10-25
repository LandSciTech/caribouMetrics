#' @include AAAClassDefinitions.R
NULL

#' Calculate binary use
#'
#' Convert continuous probability of use in each season to binary high use in
#' any season or each season using thresholds given in a table.
#' 
#'
#' @param x CaribouHabitat object where habitat use has been calculated, or a
#'   data.frame with columns PID, Spring, Summer, Fall, Winter, or a numeric
#'   vector containing ID values, in which case spring, summer, fall and winter
#'   must also be supplied
#' @param tholdTable By default \code{threshTable} which contains thresholds
#'   determimed by Rempel (2021) using Youden's J and a false
#'   negative cost of 5. Change at own risk.
#' @param bySeason logical. If FALSE (the default) the result is a single value
#'   of 1 when the habitat use was >= threshold in any season and 0 if not. If
#'   TRUE then the result is a 0 or 1 for each season.
#' @param spring,summer,fall,winter numeric vectors of habitat use probability
#'   in each season
#' @param caribouRange character. The name of the caribou range
#' @param ... passed to methods
#'
#'
#' @return
#'  If \code{bySeason} is \code{FALSE} and \code{x} is a CaribouHabitat object then the result is a RasterLayer.
#'  If \code{bySeason} is \code{TRUE} and \code{x} is a CaribouHabitat object then the result is a RasterBrick.  
#'  If \code{bySeason} is \code{FALSE} and \code{x} is a data.frame then the result is a data.frame.
#'  If \code{bySeason} is \code{TRUE} and \code{x} is a data.frame then the result is a data.frame.
#'  If \code{bySeason} is \code{FALSE} and \code{x} is a vector then the result is a vector.
#'  If \code{bySeason} is \code{TRUE} and \code{x} is a vector then the result is a data.frame.  
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
#' plot(calcBinaryUse(res))
#' plot(calcBinaryUse(res, bySeason = TRUE))
#' 
#' @export
setGeneric("calcBinaryUse", function(x, ...) standardGeneric("calcBinaryUse"))

#' @rdname calcBinaryUse
setMethod(
  "calcBinaryUse", signature(x = "CaribouHabitat"), 
  function(x, tholdTable = threshTable, bySeason = FALSE){
    caribouRange <- x@attributes$caribouRange$coefRange
    
    tTable <- threshTable %>% filter(Range == caribouRange) %>% 
      arrange(Season) %>% select(Season, Threshold)
    
    if(bySeason){
      binUse <- x@habitatUse >= tTable$Threshold
    } else {
      binUse <- any(x@habitatUse >= tTable$Threshold)
      names(binUse) <- "BinaryUse"
    }
    
    return(binUse)
  })

#' @rdname calcBinaryUse
setMethod(
  "calcBinaryUse", signature(x = "data.frame"), 
  function(x, caribouRange, tholdTable = threshTable, bySeason = FALSE){
    tTable <- threshTable %>% filter(Range == caribouRange) %>% 
      arrange(Season) %>% select(Season, Threshold)
    
    if(!all(c("PID", "Spring", "Summer", "Fall", "Winter") %in% names(x))){
      stop("names(x) must be PID, Spring, Summer, Fall, Winter")
    }
    
    out <- x %>% gather(Season, ProbUse, -PID) %>%
      left_join(tTable, by = "Season") %>% 
      mutate(binUse = ifelse(ProbUse >= Threshold, 1, 0))
    
    if(bySeason){
      out <- out %>% select(PID, Season, binUse) %>%
        pivot_wider(names_from = Season, values_from = binUse) %>% 
        rename_at(vars(-PID), ~paste0(.x, "_Cat2"))
      
      return(out)
      
    } else {
      out <- out %>% group_by(PID) %>% 
        summarise(binUse = ifelse(any(binUse == 1), 1, 0)) %>% 
        arrange(match(PID, x$PID))
    }
    
    return(out)
  })

#' @rdname calcBinaryUse
setMethod(
  "calcBinaryUse", signature(x = "numeric"), 
  function(x, spring, summer, fall, winter, caribouRange, 
           tholdTable = threshTable, bySeason = FALSE){
    
    tTable <- threshTable %>% filter(Range == caribouRange) %>% 
      arrange(Season) %>% select(Season, Threshold)
    
    x <- data.frame(PID = x, Spring = spring, Summer = summer, Fall = fall, 
                    Winter = winter)
    x <- calcBinaryUse(x, caribouRange, tholdTable, bySeason)
    
    if(bySeason){
      return(select(x, -PID))
    } else {
      
      return(x$binUse)
    }
  })