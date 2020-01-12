#' generateCovBand function, R package version of MultiTrans
#'
#' generateCovBand function estimate correlation in the rotated space
#'
#' @param windowSize is the window size defined in terms of the number of SNPs.
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
#'    gR <- generateR(X, K, VC, outPath ="./", outName = "r.txt")
#'
#'    covBand <- generateCovBand(1000, gR, outPath = "./", outName = "c.txt")
#' @export
#'
generateCovBand <- function(windowSize, gR, outPath, outName) {
  ptm <- proc.time()
  Rcpp::sourceCpp("src/generateCovBand.cpp")

  out <- paste(outPath, "/", outName, sep = "")
  generateCovBand(windowSize, gR, out)

  print(proc.time() - ptm)
  return(paste(outPath, "/", outName, sep = ""))
}
