#' generateCovBand function, R package version of MultiTrans
#'
#' generateCovBand function estimate correlation in the rotated space
#'
#' @param c is obtained from the generateCovBand function
#' @param windowSize is the window size defined in terms of the number of SNPs.
#' @param maxstat max stat file
#' @param sampling 100 used in MultiTrans paper, read MultiTrans paper for the detail, https://doi.org/10.1186/s13059-016-0903-6
#' @param seed 123 used in MultiTrans paper, read MultiTrans paper for the detail, https://doi.org/10.1186/s13059-016-0903-6
#' @param sorted storing the sampled maximum statistics over the markers.
#' @param threshold is a text file containing pointwise p-values you want to correct, delimitered by space or newline.
#' @param outPath is a parameter that specifies the path of the result file
#'
#' @return estimated correlation in the rotated space
#'
#' @examples
#'    X <- as.matrix(read.table(SNPData))
#'    K <- Kinship(X, outPath, outName)
#'    VC <- varComp(K, Y, X, outPath, outName)
#'    gR <- generateR(X, K, VC, outPath ="./", outName = "r.txt")
#'    gC <- generateCovBand(1000, gR, outPath = "./", outName = "c.txt")
#'
#'    slide <- Slide(gC, windowSize, "./maxstat", 100, 123, "./sorted", "./threshold.txt", "./result")
#' @export
#'
Slide <- function(c, windowSize, maxstat, sampling, seed, sorted, threshold, outPath){
  ptm <- proc.time()
  # Rcpp::sourceCpp("src/slide_1prep.c")

  out <- paste(outPath, "/", "prep", sep = "")

  system(paste("src/slide_1prep -C", c, windowSize, out)) # For band covariance matrix
  system(paste("src/slide_2run", out, maxstat, sampling, seed))
  system(paste("src/slide_3sort", sorted, maxstat))
  system(paste("src/slide_4correct -p", sorted, threshold, paste(outPath, "/", "MultiTrans.out", sep = "")))

  result <- as.matrix(read.table(paste(outPath, "/", "MultiTrans.out", sep = "")))
  print(proc.time() - ptm)

  return(result)
}
