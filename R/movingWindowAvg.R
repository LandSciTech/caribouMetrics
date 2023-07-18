#' Calculate the moving window average of raster
#'
#' Calculates a moving window average for each raster cell using a circular
#' window. If `veloxpkg = TRUE` it uses the velox package.
#'
#' @param rast a Raster* object which must have equal x and y resolution or they
#'   will be forced to match
#' @param radius the radius of the circular window
#' @param nms the names of each raster layer
#' @param offset Should offsetting be used to match Hornseth and Rempel 2016
#' @param na.rm,pad,padValue arguments passed on to terra::focal
#' 
#' @noRd

movingWindowAvg <- function(rast, radius, nms, offset = TRUE, 
                            naInternal = c("ignore", "interpolate", "zero"),
                            naExternal = c("ignore", "NA", "expand"),
                            usePfocal = FALSE){
  naInternal <- match.arg(naInternal)
  naExternal <- match.arg(naExternal)
  
  # changed later if needed 
  doMask <- FALSE
  mean_fun <- "mean"
  na.rm <- FALSE
  padValue <- NA
  pad <- FALSE
  doAltFun <- FALSE
  
  nl <- terra::nlyr(rast)

  cf2 <- terra::focalMat(rast, radius, "circle")
  
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
  
  if(usePfocal){
    if(!requireNamespace("pfocal", quietly = TRUE)){
      stop("Install pfocal with remotes::install_github('LandSciTech/pfocal') ",
           "to use this option")
    }
    if(!pad){
      if(na.rm){
        # ignore NAs inside the raster but still remove them on the edges
        rast <- terra::classify(rast, data.frame(NA,0), others = NULL)
        # pfocal does not have a pad argument but the equivalent is:
        padValue <- NA
        # na.rm <- FALSE
      } else {
        padValue <- NA
      }
      
    }
    if(nl == 1){
      rast <- pfocal::pfocal(rast, kernel = cf2, na.rm = na.rm, 
                             edge_value = padValue)
    } else {
      for(i in 1:nl){
        rast[[i]] <- pfocal::pfocal(rast[[i]], kernel = cf2, na.rm = na.rm, 
                                    edge_value = padValue)
      }
    }
  } else {
    if(naInternal == "interpolate"){
      rastIn <- rast
      # interpolate NAs inside the raster 
      rast <- terra::focal(rast, nrow(cf2), "mean", na.rm = TRUE, 
                           na.policy = "only")
      
      doMask <- TRUE
    } else if(naInternal == "ignore") {
      na.rm <- TRUE
      if(packageVersion("terra") < "1.7.41"){
        mean_fun <- function(i, weights, na.rm){
          weighted.mean(i, w = weights, na.rm=na.rm)
        } 
        
        w <- as.vector(cf2)
        cf2 <- ncol(cf2)
        doAltFun <- TRUE
        
        message("update to terra version > 1.7.41 for faster moving window")
      }
    } else if(naInternal == "zero"){
      rastIn <- rast
      # Set NAs inside the raster to 0
      rast <- terra::subst(rast, from = NA, to = 0)
      
      doMask <- TRUE 
    } else {
      stop("naInternal does not have an expected value", call. = FALSE)
    }
    
    if(naExternal == "NA"){
      padValue <- NA
      na.rm <- FALSE
      if(naInternal == "ignore"){
        stop("naInternal == 'ignore' and naExternal == 'NA' are not compatible",
             call. = FALSE)
      }
    } else if (naExternal == "ignore"){
      na.rm <- TRUE
    } else if( naExternal == "expand"){
      pad <- TRUE
    }
    
    if(doAltFun){
      if(terra::nlyr(rast) > 1){
        lapply(1:terra::nlyr(rast), \(x){
          rast <- terra::focal(rast[[x]], w = cf2, fun = mean_fun,
                               weights = w,
                               na.rm = na.rm, 
                               fillvalue = padValue,
                               expand = pad,
                               na.policy = "omit")
        })
      }

    } else {
      rast <- terra::focal(rast, w = cf2, fun = mean_fun,
                           na.rm = na.rm, 
                           fillvalue = padValue,
                           expand = pad,
                           na.policy = "omit")
    }

    
  }
  
  
  if(doMask){
    # add back internal NAs with masking
    rast <- terra::mask(rast, rastIn)
  }
  
  
  rast <- rast %>% `names<-`(nms)
  
  return(rast)
}

# terra focal examples for understanding inputs
   
# v <- vect(system.file("ex/lux.shp", package="terra"))
# r <- rast(system.file("ex/elev.tif", package="terra"))
# r[45:50, 45:50] <- NA
# 
# m <- focalMat(r, 0.01)
# 
# # also try "mean" or "min"
# f <- "sum" 
# # na.rm=FALSE
# plot(focal(r, m, f) , fun=lines(v))
# 
# # na.rm=TRUE
# plot(focal(r, m, f, na.rm=TRUE), fun=lines(v))
# 
# # only change cells that are NA
# plot(focal(r, m, f, na.policy="only", na.rm=TRUE), fun=lines(v))
# 
# # do not change cells that are NA
# plot(focal(r, m, f, na.policy="omit", na.rm=TRUE), fun=lines(v))
# 
# plot(focal(r, m, f, na.policy="omit", na.rm=FALSE, fillvalue = 0), fun=lines(v))
# 
# plot(focal(r, m, f, na.rm=TRUE, fillvalue = 0), fun=lines(v))
