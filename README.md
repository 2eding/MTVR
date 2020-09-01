# MTVR
R version of MultiTrans
(http://genetics.cs.ucla.edu/multiTrans/)

# Prerequistes for running Source Code
(suggest using Linux environment)

# Usage

### Step 1. Install and load package
library(devtools)<br>
install_github("2eding/MTVR")<br>
library(MTVR)<br><br>
### Step 2. Input data
X <- as.matrix(data.table::fread("X_rightdim.txt")) (individual x SNP)<br>
Y <- as.matrix(data.table::fread("Y.txt")) (Phenotype x individual)<br>
* data.table::fread <= This package is useful for input large data<br>
* X_rightdim.txt is snp data
* Y.txt is gene expression data<br><br>
### Step 3. Calculate the kinship
K <- Kinship(X, outPath, outName)<br><br>
### Step 4. Calculate the variance components
VC <- varComp(Y, K, numThreads, outPath, outName) (numThreads = Determines how many threads to use)<br><br>
### Step 5. Estimate correlation in the rotated space
corrMatrix <- generateR(X, K, VC, outPath, outName)<br>
covBandMatrix <- generateC(windowSize, corrMatrix, outPath = "./", outName = "c.txt")<br><br>
### Step 6. Run SLIDE
step1 <- system.file('./slide_1prep', package='MTVR')
step2 <- system.file('./slide_2run', package='MTVR')
step3 <- system.file('./slide_3sort', package='MTVR')
step4 <- system.file('./slide_4correct', package='MTVR')<br>
system(paste(step1, "-C", "covBandMat path", windowSize, "outPath & outName"))
system(paste(step2, "step1's outPath & outName", "outPath & outName", simulationNum, seedNum))
system(paste(step3, "outPath & outName", "step2's outPath & outName"))
system(paste(step4, "-p", "step2's outPath & outName", "threshold file path", "outPath & outName"))
