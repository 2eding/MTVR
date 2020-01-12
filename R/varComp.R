#' varComp
#'
#' varComp function will estimate variance components
#'
#' @param K is obtained from the Kinship function
#' @param Y is the phenotype matrix, individual x phenotypes
#' @param X is the SNP matrix, individual x SNPs
#' @param outPath is a parameter that specifies the path of the result file
#' @param outName is a parameter that specifies the name of the result file
#'
#' @return Variance components Vg: genetic factor Ve: environment factor
#'
#' @importFrom lmmlite eigen_rotation
#' @importFrom lmmlite fitLMM
#'
#' @examples
#'    X <- as.matrix(read.table(SNPData))
#'    Y <- as.matrix(read.table(GeneExpressionData))
#'    K <- Kinship(X, outPath, outName)
#'
#'    VC <- varComp(K, Y, X, outPath, outName)
#' @export


varComp <- function(K, Y, X, outPath, outName){
  ptm <- proc.time()
  K_origin <- K

  for (i in 1:dim(Y)[2]) {
    K <- K_origin
    print(i)

    e = lmmlite::eigen_rotation(K, t(Y)[i,], use_cpp = T)
    VC = lmmlite::fitLMM(
      e$Kva,
      e$y,
      e$X,
      reml = T,
      use_cpp = T,
      tol = 1e-6,
      check_boundary = T
    )

    write.table(
      VC$sigmasq_g,
      "Vg_temp.txt",
      row.names = F,
      col.names = F,
      append = T,
      quote = F,
      sep = "\n"
    )
    write.table(
      VC$sigmasq_e,
      "Ve_temp.txt",
      row.names = F,
      col.names = F,
      append = T,
      quote = F,
      sep = "\n"
    )
  }

  VCbind <- cbind(vg=as.matrix(read.table("Vg_temp.txt")), ve=as.matrix(read.table("Ve_temp.txt")))
  file.remove("Vg_temp.txt")
  file.remove("Ve_temp.txt")
  write.table(VCbind, paste(outPath, "/", outName, sep = ""), row.names = F, col.names = F, quote = F)
  vc <- VCbind

  print(proc.time() - ptm)
  return(vc)
}
