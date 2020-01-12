#Copyright (c) 2016 ZarLab

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' generateR function, R package version of MultiTrans
#'
#' generateR function estimate correlation in the rotated space
#'
#' @param X is the SNP matrix, individual x snp
#' @param K is obtained from the Kinship function
#' @param VC is obtained from the varComp function
#' @param outPath is a parameter that specifies the path of the result file
#' @param outName is a parameter that specifies the name of the result file
#'
#' @return estimated correlation in the rotated space
#'
#' @examples
#'    X <- as.matrix(read.table(SNPData))
#'    K <- Kinship(X, outPath, outName)
#'    VC <- varComp(K, Y, X, outPath, outName)
#'
#'    gR <- generateR(X, K, VC, outPath, outName)
#' @export

############### define functions ###############
generateR <- function(X, K, VC, outPath, outName) {
  ptm <- proc.time()
  getZscore = function(snp_one,pheno_one) {                               #for (one snp* every indi) and (one pheno * every indi)
    coeff = summary(lm(pheno_one~snp_one))$coeff
    if (dim(coeff)[1] == 1) return(NA)
    else if(is.na(coeff[2,4])) return(NA)   ### this has been changed
    else {
      zscore = abs(qnorm(coeff[2,4]/2))
      if (coeff[2,3] >= 0) return(zscore)
      else return(-1*zscore)
    }
  }
  chol_solve <- function(K) {
    a = eigen(K)$vectors
    b = eigen(K)$values
    b[b<1e-13] = 1e-13
    b = 1/sqrt(b)
    return(a%*%diag(b)%*%t(a))
  }
  rotate <- function(Y, sigma) {
    U <- chol_solve(sigma)
    tU <-t(U)
    UY = tU%*%Y
    return(UY)
  }

  snpNum <- dim(X)[2]
  indiNum <- dim(X)[1]
  Vg <- median(VC[,1])
  Ve <- median(VC[,2])
  I <- matrix(0, nrow = indiNum, ncol = indiNum)
  I[row(I) == col(I)] <- 1
  sigmaM <- Vg*K + Ve*I
  UX <- rotate(X, sigmaM)
  Ur <- cor(UX, UX)
  write.table(Ur, paste(outPath, "/", outName, sep = ""), row.names = F, col.names = F, quote = F)

  print(proc.time() - ptm)
  return(paste(outPath, "/", outName, sep = ""))
}
