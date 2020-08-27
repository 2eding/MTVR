#' varComp
#'
#' varComp function will estimate variance components
#'
#' @param Y is the phenotype matrix, phenotypes x individual
#' @param K is obtained from the Kinship function
#' @param numThreads Number of parallel processes
#' @param outPath is a parameter that specifies the path of the result file
#' @param outName is a parameter that specifies the name of the result file
#'
#' @return Variance components Vg: genetic factor Ve: environment factor
#'
#' @examples
#'    Y <- as.matrix(read.table(GeneExpressionData))
#'    K <- Kinship(X, outPath, outName)
#'
#'    VC <- varComp(Y, K, numThreads, outPath, outName)
#' @export

varComp <- function(Y, K, numThreads, outPath, outName){
  ptm <- proc.time()

  list_Y <- as.list(as.data.frame(t(Y)))
  length(list_Y)

  raw_vc <- pbmcapply::pbmclapply(list_Y, function(y){
    regress::regress(y ~ 1, ~ K, pos=c(T, T), tol=1e-8)$sigma
  },mc.cores = numThreads)

  print(proc.time() - ptm)
  VC <- do.call(rbind, raw_vc)
  VC <- as.matrix(VC)

  VC[,1] <- round(VC[,1], digits = 5)
  VC[,2] <- round(VC[,2], digits = 5)
  write.table(VC, paste(outPath, "/", outName, sep = ""), row.names = F, col.names = F, quote = F)

  return(VC)
}
