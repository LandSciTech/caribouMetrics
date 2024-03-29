% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{plcToResType}
\alias{plcToResType}
\title{Lookup table for PLC to resource type}
\format{
A data frame with 30 rows and 2 variables:
\itemize{
\item PLCCode: Provincial landcover number code
\item ResourceType: Letter code indicating the resource type
}
}
\source{
LSL script for work published in:
Hornseth, M.L. and Rempel, R.S., 2016. Seasonal resource selection of
woodland caribou (Rangifer tarandus caribou) across a gradient of
anthropogenic disturbance. Canadian Journal of Zoology, 94(2), pp.79-93.
\url{https://doi.org/10.1139/cjz-2015-0101}
}
\usage{
plcToResType
}
\description{
A table to convert Provincial Land Cover classes to resource types used in
caribou resource selection model.
}
\seealso{
Caribou habitat functions:
\code{\link{CaribouHabitat-class}},
\code{\link{calcBinaryUse}()},
\code{\link{caribouHabitat}()},
\code{\link{coefTableHR}},
\code{\link{coefTableStd}},
\code{\link{fnlcToResType}},
\code{\link{loadSpatialInputs}()},
\code{\link{rasterizeLineDensity}()},
\code{\link{reclassPLC}()},
\code{\link{resTypeCode}},
\code{\link{results}()},
\code{\link{rfuToResType}},
\code{\link{threshTable}},
\code{\link{updateCaribou}()}
}
\concept{habitat}
\keyword{datasets}
