#' Calculate the moving window average of raster
#'
#' Calculates a moving window average for each raster cell using a circular
#' window. If \code{veloxpkg = TRUE} it uses the velox package.
#'
#' @param rast a Raster* object which must have equal x and y resolution or they
#'   will be forced to match
#' @param radius the radius of the circular window
#' @param nms the names of each raster layer
#' @param veloxpkg should the velox package be used?

movingWindowAvg <- function(rast, radius, nms, offset = TRUE, 
                            na.rm = TRUE, pad = FALSE, padValue = NA){
  
  nl <- nlayers(rast)
  # if(raster::res(rast)[1] != raster::res(rast)[2]){
  #   raster::res(rast) <- c(raster::res(rast)[1], raster::res(rast)[1])
  # }
  cf2 <- focalWeight(rast, radius, "circle")
  
  if(nrow(cf2) > ncol(cf2)){
    dif <- nrow(cf2) - ncol(cf2)
    colsToAdd <- matrix(rep(0, times = nrow(cf2)*dif), nrow = nrow(cf2),
                        ncol = dif)
    cf2 <- cbind(cf2, colsToAdd)
  }
  if(ncol(cf2) > nrow(cf2)){
    dif <- ncol(cf2) - nrow(cf2)
    rowsToAdd <- matrix(rep(0, times = ncol(cf2)*dif), nrow = dif,
                        ncol = ncol(cf2))
    cf2 <- rbind(cf2, rowsToAdd)
  }
  
  if(offset > 0){
    # add rows to make only far edge points overlap for 4 horizontal/vertical
    # offsets
    nToAdd2 <- nrow(cf2)-1
    cf2_off <- matrix(0, nrow = nrow(cf2) + nToAdd2, ncol = nrow(cf2) + nToAdd2)
    
    cf2_off5 <- cf2_off
    cf2_off5[(nToAdd2/2+1):(nToAdd2/2+nrow(cf2)), (nToAdd2+1):(nToAdd2+nrow(cf2))] <- cf2
    cf2_off6 <- cf2_off
    cf2_off6[(nToAdd2+1):(nToAdd2+nrow(cf2)), (nToAdd2/2+1):(nToAdd2/2+nrow(cf2))] <- cf2
    cf2_off7 <- cf2_off
    cf2_off7[0:(nToAdd2+1), (nToAdd2/2+1):(nToAdd2/2+nrow(cf2))] <- cf2
    cf2_off8 <- cf2_off
    cf2_off8[(nToAdd2/2+1):(nToAdd2/2+nrow(cf2)), 0:(nToAdd2+1)] <- cf2
    
    # The centre "offset"
    cf2_off9 <- cf2_off
    cf2_off9[(nToAdd2/2+1):(nToAdd2/2+nrow(cf2)), (nToAdd2/2+1):(nToAdd2/2+nrow(cf2))] <- cf2
    
    # find top left edge of circle
    suppressWarnings({
      test <- FALSE
      for (i in 1:nrow(cf2)) {
        x <- which(cf2[i, ] > 0)
        y <- which(cf2[, i] > 0)
        test <- any(x == i & x == y)
        if(test){
          break
        }
      }
    })
    
    # make cell at i,i overlap point for 4 diagonal offsets
    nToAdd <- nToAdd2-i-2
    
    cf2_off1 <- cf2_off
    cf2_off1[(nToAdd+1):(nToAdd+nrow(cf2)),(nToAdd+1):(nToAdd+nrow(cf2))] <- cf2 
    cf2_off2 <- cf2_off
    cf2_off2[(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1),(nToAdd+1):(nToAdd+nrow(cf2))] <- cf2
    cf2_off3 <- cf2_off
    cf2_off3[(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1),(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1)] <- cf2
    cf2_off4 <- cf2_off
    cf2_off4[(nToAdd+1):(nToAdd+nrow(cf2)),(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1)] <- cf2
    
    if(offset == 16){
      # add 8 more diagonal offsets
      
      nToAdd3 <- nToAdd2 - (i*2 - 1)
      nToAdd4 <- nToAdd + (nToAdd - nToAdd3)/2 
      
      cf2_off10 <- cf2_off
      cf2_off10[(nToAdd3+1):(nToAdd3+nrow(cf2)),(nToAdd4+1):(nToAdd4+nrow(cf2))] <- cf2 
      cf2_off11 <- cf2_off
      cf2_off11[(nToAdd4+1):(nToAdd4+nrow(cf2)),(nToAdd3+1):(nToAdd3+nrow(cf2))] <- cf2 
      
      cf2_off14 <- cf2_off
      cf2_off14[(nrow(cf2)-nToAdd3):(nrow(cf2)*2-nToAdd3-1),(nToAdd4+1):(nToAdd4+nrow(cf2))] <- cf2 
      cf2_off15 <- cf2_off
      cf2_off15[(nrow(cf2)-nToAdd4):(nrow(cf2)*2-nToAdd4-1),(nToAdd3+1):(nToAdd3+nrow(cf2))] <- cf2 
      
      cf2_off12 <- cf2_off
      cf2_off12[(nrow(cf2)-nToAdd4):(nrow(cf2)*2-nToAdd4-1),(nrow(cf2)-nToAdd3):(nrow(cf2)*2-nToAdd3-1)] <- cf2
      cf2_off13 <- cf2_off
      cf2_off13[(nToAdd3+1):(nToAdd3+nrow(cf2)),(nrow(cf2)-nToAdd4):(nrow(cf2)*2-nToAdd4-1)] <- cf2
      
      cf2_off16 <- cf2_off
      cf2_off16[(nrow(cf2)-nToAdd3):(nrow(cf2)*2-nToAdd3-1),(nrow(cf2)-nToAdd4):(nrow(cf2)*2-nToAdd4-1)] <- cf2
      cf2_off17 <- cf2_off
      cf2_off17[(nToAdd4+1):(nToAdd4+nrow(cf2)),(nrow(cf2)-nToAdd3):(nrow(cf2)*2-nToAdd3-1)] <- cf2
      
      cf2 <- cf2_off1 + cf2_off2 + cf2_off3 + cf2_off4 + cf2_off5 + cf2_off6 + 
        cf2_off7 + cf2_off8 + cf2_off10 + cf2_off11+ cf2_off12 + cf2_off13 +
        cf2_off14 + cf2_off15 + cf2_off16 + cf2_off17
      cf2 <- cf2/16
    } else {
      cf2 <- cf2_off1 + cf2_off2 + cf2_off3 + cf2_off4 + cf2_off5 + cf2_off6 + 
        cf2_off7 + cf2_off8 + cf2_off9
      cf2 <- cf2/9
    }
    
  }
  
  # error if rast is smaller than the focal window
  if(any(!dim(rast)[1:2]>dim(cf2))){
    stop("The project area is smaller than the window area", call. = FALSE)
  }
  
  if(nl == 1){
    rast <- raster::focal(rast, w = cf2, na.rm = na.rm, pad = pad,
                          padValue = padValue)
  } else {
    for(i in 1:nl){
      rast[[i]] <- raster::focal(rast[[i]], w = cf2, na.rm = na.rm, pad = pad,
                                 padValue = padValue)
    }
  }
  
  rast <- rast %>% `names<-`(nms)
  
  return(rast)
}

   