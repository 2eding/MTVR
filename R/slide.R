#' slide
#'
#' slide function will execute SLIDE
#' see detail (http://slide.cs.ucla.edu/index.html)
#'
#' @param covBandMatrix is the SNP matrix, individual x SNPs
#' @param windowSize is the window size defined in terms of the number of SNPs.
#' @param numSampling is the number of sampling.
#' @param seedNum is set seed.
#' @param threshold is a text file containing pointwise p-values you want to correct, delimitered by space or newline.
#' @param outPath is a parameter that specifies the path of the result file
#' @param outName is a parameter that specifies the name of the result file
#'
#' @return kinship coefficient
#'
#' @examples
#'    X <- as.matrix(read.table(SNPData))
#'    K <- Kinship(X, outPath, outName)
#'    VC <- varComp(Y, K, numThreads, outPath, outName)
#'    corrMatrix <- generateR(X, K, VC, outPath, outName)
#'    covBandMatrix <- generateC(windowSize, corrMatrix, outPath = "./", outName = "c.txt")
#'
#'    result <- slide(covBandMatrix, windowSize, numSampling, seedNum, threshold, outPath, outName)
#' @export
#'
# 5. run slide (You can download slide from http://slide.cs.ucla.edu/)
# 5.1. slide_1prep
# ./slide.1.0/slide_1prep -C ./testData/c.txt 1000 ./testData/prep
# 5.2. slide_2run
# ./slide.1.0/slide_2run ./testData/prep ./testData/maxstat 100 123
# 5.3. slide_3sort
# ./slide.1.0/slide_3sort ./testData/sorted ./testData/maxstat
# 5.4. slide_4correct
# ./slide.1.0/slide_4correct -p ./testData/sorted threshold.txt ./testData/MultiTrans.output

slide <- function(covBandMatrix, windowSize, numSampling, seedNum, threshold, outPath, outName){
  slide_step1 <- system.file('slide/slide_1prep', package='MTVR')
  Sys.chmod(slide_step1)
  system(slide_step1)

  slide_step2 <- function(){
    source(system.file('slide/slide_2run', package='MTVR'))
  }
  slide_step3 <- function(){
    source(system.file('slide/slide_3sort', package='MTVR'))
  }
  slide_step4 <- function(){
    source(system.file('slide/slide_4correct', package='MTVR'))
  }

  slide_step1(paste("-C", covBandMatrix, windowSize, outPath, sep = ""))

  return(result)
}
