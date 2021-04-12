###############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Script name: Real Data Demonstration
#
# Purpose of script: This script implements a simple example of the 
#                    caribouHabitat() function from the caribouMetrics package.
#                    The script uses empirical data collected for the Ring of
#                    Fire area.
#
# Author: Dr. Craig Simpkins (simpkinscraig063@gmail.com)
#
# Date created: 09/04/2021
#
# Notes: This script fundamentally is a copy of the demoMultRanges.R script
#        written by Sarah Endicott to demonstrate how to run caribouHabitat()
#        for multiple ranges simultaneously. At the time of writing the 
#        functionality for multiple ranges is only in the issue-2-multirange
#        git branch of the package. Additionally, only the data required for
#        this example is used here and not all available data in the Ring of
#        Fire directory
#
#        Lastly, note that running this script, particularly on multiple ranges
#        can be quite slow taking 15+ minutes, for that reason the multiple
#        ranges call is commented out.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############################################################################

#~~SETUP~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set the path to the data (this will need to be changed for different users)
pth_base <- "../Development_Scripts/Data_for_RSPF_models/inputNV/"

devtools::load_all(".") # Or library(caribouMetrics)

# Create a list of all the range names in the RoF.
caribouRanges <- c("Pagwachuan", "Missisa", "Ozhiski", "Nipigon", "James Bay")
caribouRangeCoefs <- rev(caribouRanges)

#~~INPUT DATA~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the shape file containing the boundary for all the caribou ranges
# stored in the RoF directory and then filter down to the named ranges
projectPoly <- st_read(paste0(pth_base,
                              "caribouRanges/Caribou_Range_Boundary.shp")) %>% 
  filter(RANGE_NAME %in% caribouRanges) %>% 
  rename(Range = RANGE_NAME)

# Import the provincial landcover data and reclassify its values to resource 
# types. Then lower the resolution down to that matching the original
# Hornseth and Rempel model.
landCover <- raster(paste0(pth_base, "ROFData/plc50.tif")) %>% 
  reclassPLC() %>% 
  raster::aggregate(5, fun = raster::modal)

esker <- st_read(paste0(pth_base, "ROFData/esker.shp"))

# Import linear features and store them as a list object for use in the
# caribouHabitat() function. 
linFeat <- list(roads = st_read(paste0(pth_base, 
                                       "ROFData/road_ORNMNRFROF2010.shp")),
                rail = st_read(paste0(pth_base, "ROFData/rail.shp")),
                utilities = st_read(paste0(pth_base, "ROFData/util2010.shp")))

# Create a two column data frame to allow users to specify the ranges and the
# coefficients for each range if multiple ranges are run simultanuously
caribouRange <- data.frame(Range = caribouRanges, 
                           coefRange = caribouRanges, stringsAsFactors = FALSE)


#~~RUN FUNCTION~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## For multiple ranges

# I have commented this out at this stage to prevent excessive run times.
# sameCoefMultRange <- caribouHabitat(landCover = landCover,
#                                     esker = esker, 
#                                     linFeat = linFeat,  
#                                     projectPoly = projectPoly, 
#                                     caribouRange = caribouRange, 
#                                     padProjPoly = TRUE)
# 
# 
# plot(sameCoefMultRange)

## For single range (as has been standard up to this point)
SRange <- caribouHabitat(landCover = landCover,
                         esker = esker, 
                         linFeat = linFeat,  
                         projectPoly = projectPoly[projectPoly$Range 
                                                   == "Nipigon",], 
                         caribouRange = "Nipigon", 
                         padProjPoly = TRUE)

plot(SRange)
