
Slide <- function(){
  ptm <- proc.time()
  Rcpp::sourceCpp("src/generateCovBand.cpp")

  out <- paste(outPath, "/", outName, sep = "")
  generateCovBand(windowSize, gR, out)

  covBand <- as.matrix(read.table(paste(outPath, "/", outName, sep = "")))
  print(proc.time() - ptm)
  return()
}
