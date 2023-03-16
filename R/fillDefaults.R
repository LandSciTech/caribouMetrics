fillDefaults <- function(scns = NULL,
                         iF = 0, iA = 0, aS = 0, aSf = 4,
                         rS = 1, sS = 1,
                         rQ = 0.5, sQ = 0.5, J = 20, P = 1, N0 = 1000,
                         adjustR = F, assessmentYrs = 1,
                         ri = NA, cmult = NA, cw = NA,
                         curYear = 2023) {
  defList <- c(as.list(environment()))
  defList$scns <- NULL
  if (is.null(scns)) {
    scns <- as.data.frame(defList)
  } else {
    # keep all values in scns and add any that are missing using values in
    # defList
    scns <- cbind(scns, defList[which(!names(defList) %in% names(scns))])
  }

  # remove columns that are all NA because they should be missing and order like
  # defList
  scns <- select(scns, all_of(names(defList)), -where(~all(is.na(.x))))

  if (is.element("cmult", names(scns)) & is.element("cw", names(scns))) {
    stop("Specify number of cows per year in recruitment survey (cw) or multiplier of number of collared cows in recruitment survey (cmult), but not both.")
  }
  scns$ID <- seq(1:nrow(scns))
  scns$label <- ""
  for (n in names(scns)[(length(names(scns)) - 1):1]) {
    scns$label <- paste0(scns$label, n, scns[[n]], "_")
  }

  if (!is.element("iYr", names(scns))) {
    scns$iYr <- scns$curYear - scns$P + 1
  }
  return(scns)
}
