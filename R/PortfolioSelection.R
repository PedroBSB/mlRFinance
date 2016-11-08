#' Execute the SVM model with cross-validation.
#' Binary Selection
#' @param train A training dataframe with the binary values (-1 or +1)
#' @param valid A validation dataframe with the returns.
#' @param C Vector of costs to be validated.
#' @param kernel Type of kernel to be used.
#' @param parmMat Matrix of parameters of the kernel.
#' @param typePredict 0-Binary(-1,+1), 1-Probability, 2- Raw result
#' @return Validation results of the SVM for Portfolio Selection
#' @examples
#' add(1, 1)
#' add(10, 1)
PortfolioSVML1 <- function(train, valid, C, kernel, parmMat, typePredict) {
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("foreach needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("doParallel needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("parallel needed for this function to work. Please install it.",
         call. = FALSE)
  }

  #Combination with parms
  matAll<-merge(C,parmMat)

  #Insert names
  colnames(matAll)<-c("C", paste0("Parm",seq(0,ncol(parmMat)-1)))

  #Get the total number of cores
  ncl<- parallel::detectCores()

  #Register the clusters
  cl <- parallel::makeCluster(ncl)
  doParallel::registerDoParallel(cl)

  #Initialize the validation
  svmPort <- foreach(i=1:nrow(matAll), .combine=rbind, .errorhandling='pass', .packages="mlRFinance") %dopar% {
    #Cost
    C0<-matAll$C[i]
    #Parms Mean
    parmsM<-as.numeric(matAll[i,2:ncol(matAll)])
    #Training the machine
    svm<-PortfolioSelectionCSVML1(train, valid, C0,kernel, parmsM,typePredict )
    if(typePredict==0){
      res<-data.frame(matAll[i,],
                      "Prevalence" = svm$ErrorMeasureValidation$Prevalence,
                      "TruePositiveRate" = svm$ErrorMeasureValidation$TruePositiveRate,
                      "FalseNegativeRate" = svm$ErrorMeasureValidation$FalseNegativeRate,
                      "FalsePositiveRate" = svm$ErrorMeasureValidation$FalsePositiveRate,
                      "TrueNegativeRate" = svm$ErrorMeasureValidation$TrueNegativeRate,
                      "Accuracy" = svm$ErrorMeasureValidation$Accuracy,
                      "PositivePredictiveValue" = svm$ErrorMeasureValidation$PositivePredictiveValue,
                      "FalseOmissionRate" = svm$ErrorMeasureValidation$FalseOmissionRate,
                      "FalseDiscoveryRate" = svm$ErrorMeasureValidation$FalseDiscoveryRate,
                      "NegativePredictiveValue" = svm$ErrorMeasureValidation$NegativePredictiveValue,
                      "PositiveLikelihoodRatio" = svm$ErrorMeasureValidation$PositiveLikelihoodRatio,
                      "NegativeLikelihoodRatio" = svm$ErrorMeasureValidation$NegativeLikelihoodRatio,
                      "DiagnosticOddsRatio" = svm$ErrorMeasureValidation$DiagnosticOddsRatio)

    }
    else if(typePredict==1){
      res<-data.frame(matAll[i,],"MSE"=svm$ErrorMeasureValidation$MSE)
    }
    else{

    }

    res
  }

  #Stop clusters
  stopCluster(cl)
  return(svmPort)
}
