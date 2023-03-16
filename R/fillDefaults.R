fillDefaults <- function(scns = NULL,
                         defList = list(
                           iF = 0, iA = 0, aS = 0, aSf = 4,
                           rS = 1, sS = 1,
                           rQ = 0.5, sQ = 0.5, J = 20, P = 1, N0 = 1000,
                           adjustR = F, ri = NULL, assessmentYrs = 1
                         ), curYear = 2023) {
  if (is.null(scns)) {
    scns <- as.data.frame(defList)
  } else {
    fillSet <- setdiff(names(defList), names(scns))

    for (i in fillSet) {
      scns[[i]] <- defList[[i]]
    }
  }
  if (is.element("cmult", names(scns)) & is.element("cw", names(scns))) {
    stop("Specify number of cows per year in recruitment survey (cw) or multiplier of number of collared cows in recruitment survey (cmult), but not both.")
  }
  scns$ID <- seq(1:nrow(scns))
  scns$label <- ""
  for (n in names(scns)[(length(names(scns)) - 1):1]) {
    scns$label <- paste0(scns$label, n, scns[[n]], "_")
  }

  if (!is.element("iYr", names(scns))) {
    scns$iYr <- curYear - scns$P + 1
  }
  return(scns)
}
