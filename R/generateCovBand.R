#' generateCovBand function, R package version of MultiTrans
#'
#' generateCovBand function estimate correlation in the rotated space
#'
#' @param windowSize 1000 used in MultiTrans paper, read MultiTrans paper for the detail, https://doi.org/10.1186/s13059-016-0903-6
#' @param gR is obtained from the generateR function
#' @param outPath is a parameter that specifies the path of the result file
#' @param outName is a parameter that specifies the name of the result file
#'
#' @return estimated correlation in the rotated space
#'
#' @examples
#'    X <- as.matrix(read.table(SNPData))
#'    K <- Kinship(X, outPath, outName)
#'    VC <- varComp(K, Y, X, outPath, outName)
#'    r <- generateR(X, K, VC, outPath ="./", outName = "r.txt")
#'
#'    covBand <- generateCovBand(1000, "r.txt", outPath = "./", outName = "c.txt")
#' @export
#'
generateCovBand <- function(windowSize, gR, outPath, outName) {
  ptm <- proc.time()
  Rcpp::sourceCpp("src/generateCovBand.cpp")

  out <- paste(outPath, "/", outName, sep = "")
  generateCovBand(windowSize, gR, out)

  covBand <- as.matrix(read.table(paste(outPath, "/", outName, sep = "")))
  print(proc.time() - ptm)
  return(covBand)
}
