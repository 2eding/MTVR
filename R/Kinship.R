#' Kinship
#'
#' Kinship function will estimate kinship coefficient
#'
#' @param X is the SNP matrix, individual x SNPs
#' @param outPath is a parameter that specifies the path of the result file
#' @param outName is a parameter that specifies the name of the result file
#'
#' @return kinship coefficient
#'
#' @examples
#'    X <- as.matrix(read.table(SNPData))
#'    K <- Kinship(X, outPath, outName)
#' @export
Kinship <- function(X, outPath, outName){
  ptm <- proc.time()

  sd <- function(data){
    return(sum((mean(data) - data)^2)/length(data))
  }

  W <- X
  n <- nrow(W)
  m <- ncol(W)
  keep <- c()

  for (i in 1:m){
    mn <- mean(W[!is.na(W[,i]),i])
    W[is.na(W[,i]),i] <- mn
    vr <- sd(W[,i])

    if (vr == 0) s <- 1.0
    else s <- sqrt(vr)

    keep <- c(keep,i)
    W[,i] <- (W[,i] - mn) / s
  }
  W <- W[, keep]
  K <- tcrossprod(W) * 1.0/m

  write.table(K, paste(outPath, "/", outName, sep = ""), row.names = F, col.names = F, quote = F)
  print(proc.time() - ptm)
  return(K)
}
