#' generateC
#'
#' generateC function will generate CovBand matrix
#'
#' @param windowSize is the window size defined in terms of the number of SNPs
#' @param corrMatrix is obtained from the generateR function
#' @param numThreads Number of parallel processes
#' @param outPath is a parameter that specifies the path of the result file
#' @param outName is a parameter that specifies the name of the result file
#'
#' @return estimated correlation in the rotated space
#'
#' @examples
#'    X <- as.matrix(read.table(SNPData))
#'    K <- Kinship(X, outPath, outName)
#'    VC <- varComp(Y, K, numThreads, outPath, outName)
#'    corrMatrix <- generateR(X, K, VC, outPath, outName)
#'
#'    covBand <- generateC(windowSize, corrMatrix, outPath = "./", outName = "c.txt")
#' @export
generateC <- function(windowSize, corrMatrix, outPath, outName){
  ptm <- proc.time()

  ret_mat <- matrix(0, nrow = nrow(corrMatrix)-1, ncol = windowSize)

  for(i in 1:(nrow(ret_mat))){
    for(j in 1:i){
      ret_mat[i, ncol(ret_mat) + j - i] <- corrMatrix[i+1, j]
    }
  }

  write.table(ret_mat, paste(outPath, "/", outName, sep = ""), row.names = F, col.names = F, quote = F)

  print(proc.time() - ptm)
  return(ret_mat)
}
