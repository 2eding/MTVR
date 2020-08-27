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
* Y.txt is gene expression data

### Step 3. Calculate the kinship
K <- Kinship(X, outPath, outName)<br><br>

### Step 4. Calculate the variance components
VC <- varComp(Y, K, numThreads, outPath, outName) (numThreads = Determines how many threads to use)<br><br>

### Step 5. Estimate correlation in the rotated space
corrMatrix <- generateR(X, K, VC, outPath, outName)<br>
covBandMatrix <- generateC(windowSize, corrMatrix, outPath = "./", outName = "c.txt")<br>

### Step 6. Run SLIDE
#slide_1prep: data pre-processing.<br>
slideone <- slide_1prep
#slide_2run: run the actual sampling.<br>
slidetwo <- slide_2run
#slide_3sort: sort the maximum statistic.<br>
slidethree <- slide_3sort
#slide_4correct: correct p-values.<br>
slidefour <- slide_4correct
