#' Execute the Garch SVR model with cross-validation.
#'
#' @param train A training dataframe with the returns.
#' @param valid A validation dataframe with the returns.
#' @param Cm Vector of costs to be validated. Mean Equation.
#' @param epsilonM Vector of Insentitive band.  Mean Equation.
#' @param kernel Type of kernel to be used - Mean Equation.
#' @param parmMat Matrix of parameters of the kernel - Mean Equation.
#' @param Cg Vector of costs to be validated. Mean Equation.
#' @param eg Vector of Insentitive band.  Mean Equation.
#' @param kernelGarch Type of kernel to be used - Garch Equation.
#' @param parmMatGarch Matrix of parameters of the kernel - Garch Equation.
#' @return Validation results of the Garch SVR
#' @examples
#' add(1, 1)
#' add(10, 1)
GarchSVR <- function(train,valid,Cm,epsilonM,kernel,parmMat,Cg,epsilonG,kernelGarch,parmMatGarch) {
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

  #Create the matrix with all combinations between C, epsilon e parms
  matAll <- expand.grid(Cm,epsilonM,Cg,epsilonG)

  #Combination with parms (Mean Equation)
  matAll<-merge(matAll,parmMat)

  #Insert names
  colnames(matAll)<-c("Cm","epsilonM","Cg","epsilonG",
                      paste0("ParmM",seq(0,ncol(parmMat)-1)))

  #Parameters
  matAll <- merge(matAll,parmMatGarch)

  #Insert names
  colnames(matAll)<-c("Cm","epsilonM","Cg","epsilonG",
                      paste0("ParmM",seq(0,ncol(parmMat)-1)),
                      paste0("ParmG",seq(0,ncol(parmMatGarch)-1)))

  #Get the total number of cores
  ncl<- parallel::detectCores()

  #Register the clusters
  cl <- parallel::makeCluster(ncl)
  doParallel::registerDoParallel(cl)

  #Initialize the validation
  svrGarch <- foreach(i=nrow(matAll), .combine=rbind, .errorhandling='pass', .packages="mlRFinance") %dopar% {
    #Cost - Mean Equation
    C0m<-matAll$Cm[i]
    #Epsilon - Mean Equation
    eps0m<-matAll$epsilonM[i]
    #Parms Mean
    parmsM<-as.numeric(matAll[i,5:(4+ncol(parmMat))])
    #Cost - Garch Equation
    C0g<-matAll$Cg[i]
    #Epsilon - Garch Equation
    eps0g<-matAll$epsilonG[i]
    #Parms Garch
    parmsG<-as.numeric(matAll[i,(5+ncol(parmMat)):(4+ncol(parmMat)+ncol(parmMat))])
    #Training the machine
    svr<-GARCHCSVRL1(train, valid, C0m, eps0m,
                     C0g, eps0g,
                     kernel, parmsM,
                     kernelGarch, parmsG)

    res<-data.frame(matAll[i,],"MSEm"=svm$ErrorMeasureValidation,"MSEg"=svm$ErrorMeasureValidationGarch)
  }

  #Stop clusters
  stopCluster(cl)
}


