% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fnlcToResType}
\alias{fnlcToResType}
\title{Lookup table for FNLC to resource type}
\format{
A data frame with 26 rows and 2 variables:
\itemize{
\item PLCCode: FNLC number code
\item ResourceType: Letter code indicating the resource type
}
}
\usage{
fnlcToResType
}
\description{
A table to convert Far North Land Cover classes to resource types used in
caribou resource selection model.
}
\details{
The FNLC classes were linked to resource types by comparing the descriptions
of FNLC classes to the PLC classes and then linking them to the corresponding
resource types
}
\seealso{
Caribou habitat functions:
\code{\link{CaribouHabitat-class}},
\code{\link{calcBinaryUse}()},
\code{\link{caribouHabitat}()},
\code{\link{coefTableHR}},
\code{\link{coefTableStd}},
\code{\link{loadSpatialInputs}()},
\code{\link{plcToResType}},
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
