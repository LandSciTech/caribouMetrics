# Compare results with different disturbance datasets
devtools::load_all()
library(tmap)
pthBase <- "inputNV/Disturbance/"

# BEAD #========================================================================
# more info available at:
# https://open.canada.ca/data/en/dataset/a71ab99c-6756-4e56-9d2e-2a63246a5e94

pthBEAD2010 <- paste0(pthBase, "BEADwithFire/Churchill_2010/Churchill_30m.gdb")

st_layers(pthBEAD2010)

# The 500 refers to the polygons and lines with a 500m buffer. The withFire
# means fire polygons based on NRCAN National Large Fire Database have been
# merged with anthro disturbance. For both all types of disturbance are
# dissolved together. So I used the original polygon data. 

beadPoly2010 <- st_read(pthBEAD2010, layer = "Churchill_30m_Disturb_Perturb_Poly")
beadLine2010 <- st_read(pthBEAD2010, layer = "Churchill_30m_Disturb_Perturb_Line")

beadPoly2010 <- beadPoly2010 %>% st_make_valid() 

# Rasterize polygons with plc as target. Don't need linear features b/c they are
# in linFeat data set
plcD <- raster("inputNV/intermediaryData/plc_aligned.tif")

beadPoly2010Rast <- raster::rasterize(beadPoly2010 %>% 
                                             st_transform(crs = st_crs(plcD)), 
                                           plcD, 
                                           field = "Class",
                                           background = 0)

# class look up for raster
beadPoly2010ClassLU <- data.frame(Class = beadPoly2010$Class %>% levels()) %>% 
  mutate(Code = 1:n())

# Bead fire data all fires from 1970 to 2010 based on NRCAN national large fire
# database
beadFire2010 <- st_read(pthBEAD2010, layer = "Churchill_Fire_40yr_1970_2010")

beadFire2010Rast <- fasterize::fasterize(beadFire2010 %>% 
                                           st_transform(crs = st_crs(plcD)), 
                                         plcD, 
                                         field = NULL,
                                         background = 0)
raster::writeRaster(beadFire2010Rast, "inputNV/intermediaryData/beadFire2010.tif",
                    overwrite = TRUE)

# save beadPolyRast as binary 
raster::writeRaster(beadPoly2010Rast > 0,
                    "inputNV/intermediaryData/beadPoly2010.tif", overwrite = TRUE)


# NRCAN National Large Fire Database #==========================================
# Downloaded National Fire Database fire polygon data (large fires only ≥ 200
# hectares) from:
# https://cwfis.cfs.nrcan.gc.ca/datamart/download/nfdbpoly?token=ab402faa420fd76a27c5c9f3c9fa9aff

# ** Not currently needed but could be useful in future **

# lgFiresNFDB <- st_read(paste0(pthBase,
#                           "NFDB_poly_large_fires/NFDB_poly_20200711_large_fires.shp"))
# 
# projectPolyD <- st_read("./inputNV/caribouRanges/Caribou_Range_Boundary.shp",
#                         quiet = TRUE) %>%
#   filter(RANGE_NAME == "Churchill")
# 
# CHlgFiresNFDB <- lgFiresNFDB %>% st_transform(crs = st_crs(projectPolyD)) %>%
#   st_make_valid() %>%
#   st_intersection(projectPolyD)
# 
# 
# st_write(CHlgFiresNFDB,
#          paste0(pthBase, "NFDB_poly_large_fires/NFDB_large_fires_Churchill.shp"),
#          append = FALSE)
# 
# rm(lgFiresNFDB)
# 
# CHlgFiresNFDB <- st_read(paste0(pthBase,
#                                 "NFDB_poly_large_fires/NFDB_large_fires_Churchill.shp"))
# 
# NFDBFire2011 <- fasterize::fasterize(CHlgFiresNFDB %>% filter(YEAR == 2011) %>% 
#                                            st_transform(crs = st_crs(plcD)), 
#                                          plcD, 
#                                          field = NULL,
#                                          background = 0)
# raster::writeRaster(beadFire2010Rast + NFDBFire2011 > 0, 
#                     "inputNV/intermediaryData/beadFire2010and2011NFDB.tif",
#                     overwrite = TRUE)
# tmap_mode("plot")
# qtm(CHlgFiresNFDB2010)
# qtm(plcD == 8)
# 
# # Fires from 1990 to 2000 correspond closely with plc depletion - burn class
# 
# # According to BEAD docs fire data from 1975 to 2011 was merged with
# # Landsat interpreted Anthro disturbance
# CHlgFiresNFDB2010 <- CHlgFiresNFDB %>% filter(YEAR < 2011, YEAR >= 1970)
# 
# tmap_mode("view")
# qtm(beadFire2010)+
#   qtm(CHlgFiresNFDB2010, fill = "YEAR")
# 
# # CHlgFiresNFDB2010 lines up well with beadFire2010 as expected could use NFDB
# # to get same data for different set of years. There are some small fires in
# # beadFire2010 that are not in NFDB but might be if we downloaded the whole
# # thing instead of just large fires.

# OFRI Disturbance database #===================================================
pthBase2 <- paste0(pthBase, "OFRIDisturbanceHistory/ChurchillCaribouBorealDisturbance")

pths <- list.files(pthBase2, pattern = ".shp$", full.names = TRUE)

plcD <- raster("inputNV/intermediaryData/plc_aligned.tif")

# binding them together messud up CRS so read one to get correct
OFRIDist1972 <- pths[1] %>% st_read()

# load all the points into one df
OFRIDistHist <- map(pths, st_read) %>% rbind_list()%>% 
  st_as_sf() %>% 
  st_set_crs(st_crs(OFRIDist1972)) %>% 
  st_transform(st_crs(plcD))

# names have been shortened re name based on documentation
colnames(OFRIDistHist) <- c("OBJECTID", "pointid", "Disturbed","DistYear", 
                            "FirLandsat", "HarLandsat", "FirAFFES", "HarMNRF", 
                            "FirMNRFSRB", "HarMNRFSRB", "FirNRCAN", "HarNRCAN",
                            "FirPrev", "HarPrev", "DistCodeLS", "FirConf", 
                            "HarConf", "geometry")

OFRIDistHist$geometry %>% st_is_empty() %>% sum()
# There are 176226 points with empty geometries. I don't know why or what the
# means. Could it be realted to subsetting the Churchill data?

# There are 4 sets of Fire and Harvest to pull out of here and rasterize
# 1) "FirLandsat", "HarLandsat" 2) "FirAFFES", "HarMNRF", 
# 3) "FirMNRFSRB", "HarMNRFSRB" 4) "FirNRCAN", "HarNRCAN"

dataSetNms <- list(OFRID = c("FirLandsat", "HarLandsat"), 
                   MNRFAFFES = c("FirAFFES", "HarMNRF"), 
                   MNRFSRB = c("FirMNRFSRB", "HarMNRFSRB"), 
                   NRCAN = c("FirNRCAN", "HarNRCAN"))


OFRIDistHist <- OFRIDistHist %>% 
  filter(!st_is_empty(geometry)) %>%
  mutate(geometry = st_cast(geometry, "POINT"))

# what years to include and what to do about points with dist in more than one year
# should be 1990 to 2010 because plc uses dist in last 10 years.
# Keep only the most recent disturbance but do it separately for each source
OFRIDistHist90_10 <- OFRIDistHist %>% filter(between(DistYear, 1990, 2010)) 

# check resolution vs points
# plot(plcD)
# ext <- raster::drawExtent()
# 
# testOFRID <- st_crop(OFRIDistHist90_10, ext)
# testplc <- raster::crop(plcD, ext)

# the points don't align perfectly with plc grid. could convert to raster of
# more similar grid first and then resample but resampling wouldn't be very
# different from current

rastizeDataSets <- function(df, nms, tarRast, thold){
  dfFire <- select(df, all_of(nms[1]), DistYear, pointid) %>% 
    filter_at(vars(1), any_vars(. > 0)) 
  
  if(nrow(dfFire) == 0){
    fire <- tarRast %>% raster::setValues(0)
    
  } else {
    dfFire <- dfFire %>%  
      arrange(desc(DistYear)) %>% # distinct keeps first row
      distinct(pointid, .keep_all = TRUE) 
    
    fire <- raster::rasterize(dfFire, tarRast, field = nms[1], background = 0, 
                              fun = "count")
    
    fire <- fire > thold
  }
  
  dfHarv <- select(df, all_of(nms[2]), DistYear, pointid) %>% 
    filter_at(vars(1), any_vars(. > 0)) 
  
  if(nrow(dfHarv) == 0){
    harv <- tarRast %>% raster::setValues(0)
    
  } else {
    dfHarv <- dfHarv %>%  
      arrange(desc(DistYear)) %>% # distinct keeps first row
      distinct(pointid, .keep_all = TRUE)
    
    harv <- raster::rasterize(dfHarv, tarRast, field = nms[2], background = 0,
                              fun = "count")
    harv <- harv > thold
  }
  
  return(lst(fire, harv))
}

# Used count to determine the number of points in the raster cell which each
# represents 1.44 ha but we need a binary disturbance variable so set threshold
# to > 2 points because that is about half 1.44*10000*2/250^2
dataSets <- purrr::map(dataSetNms, 
                       ~rastizeDataSets(OFRIDistHist90_10, .x, plcD, 2))
beepr::beep()
# write them all to files
pthWrite <- "inputNV/intermediaryData/"

wrRastFun <- function(rast, setNm, layerNm, year){
 foldNm <- paste0(pthWrite, setNm)
 
 fileNm <- paste0(foldNm, "_", layerNm, year, ".tif")
 raster::writeRaster(rast, fileNm, overwrite = TRUE)
}

purrr::walk2(1:length(dataSets) %>% rep(2), 
             c(rep(1,length(dataSets)), rep(2, length(dataSets))), 
             ~wrRastFun(rast = dataSets[[.x]][[.y]], names(dataSets)[.x], 
                        names(dataSets[[.x]])[.y], 2010))

# repeat for 2020
OFRIDistHist00_20 <- OFRIDistHist %>% filter(between(DistYear, 2000, 2020)) 

dataSets00_20 <- purrr::map(dataSetNms, 
                            ~rastizeDataSets(OFRIDistHist00_20, .x, plcD, 2))
beepr::beep()

purrr::walk2(1:length(dataSets00_20) %>% rep(2), 
             c(rep(1,length(dataSets00_20)), rep(2, length(dataSets00_20))), 
             ~wrRastFun(rast = dataSets00_20[[.x]][[.y]], names(dataSets00_20)[.x], 
                        names(dataSets00_20[[.x]])[.y], 2020))
