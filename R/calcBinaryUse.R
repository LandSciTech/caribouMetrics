#' @include AAAClassDefinitions.R
NULL

#' Calculate binary use
#'
#' Convert continuous probability in each season of use to binary high use in
#' any season using thresholds given in a table.
#'
#' @param x CaribouHabitat object where habitat use has been calculated, or a 
#' data.frame with columns PID, Spring, Summer, Fall, Winter
#' @param caribouRange The caribou range of the CaribouHabitat object
#' @param tholdTable By default \code{threshTable} which contains thresholds 
#'   determimed by Rempel and Hornseth (2018) using Youden's J and a false 
#'   negative cost of 5. Change at own risk.
#'
#' @return a binary RasterLayer with 1 for habitat that is above the 
#'   threshold in any season 
#'   
#' @export
#'
#' @examples
setGeneric("calcBinaryUse", function(x, ...) standardGeneric("calcBinaryUse"))

#' @rdname calcBinaryUse
setMethod(
  "calcBinaryUse", signature(x = "CaribouHabitat"), 
  function(x, caribouRange, tholdTable = threshTable, bySeason = FALSE){
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