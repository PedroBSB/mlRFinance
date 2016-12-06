#' Execute the Single Causality SVM (SiCSVM)
#'
#' @param ytrain A training dataframe with the dependent variable
#' @param Xtrain A training dataframe with the covariates.
#' @param wtrain A training dataframe with the binary group (1=Treated, 0=Control)
#' @param yvalid A validation dataframe with the dependent variable
#' @param Xvalid A validation dataframe with the covariates.
#' @param wvalid A validation dataframe with the binary group (1=Treated, 0=Control)
#' @param C Vector of costs to be validated.
#' @param epsilon Vector of Insentitive band.
#' @param kernel Type of kernel to be used.
#' @param parmMat Matrix of parameters of the kernel.
#' @return Causality evaluation.
#' @examples
#' add(1, 1)
#' add(10, 1)
SiCSVM <- function(ytrain,Xtrain,wtrain,yvalid,Xvalid,wvalid,C,epsilon,kernel,parmMat) {
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
  matAll <- expand.grid(C,epsilon)

  #Combination with parms (Mean Equation)
  matAll<-merge(matAll,parmMat)

  #Insert names
  colnames(matAll)<-c("C","epsilon",
                      paste0("ParmM",seq(0,ncol(parmMat)-1)))


  #Insert names
  colnames(matAll)<-c("Cm","epsilonM",
                      paste0("ParmM",seq(0,ncol(parmMat)-1)))

  #Get the total number of cores
  ncl<- parallel::detectCores()

  #Register the clusters
  cl <- parallel::makeCluster(ncl)
  doParallel::registerDoParallel(cl)

  #Initialize the validation
  svrCausality <- foreach(i=1:nrow(matAll), .combine=rbind, .errorhandling='pass', .packages="mlRFinance") %dopar% {
    #Cost - Mean Equation
    C0m<-matAll$C[i]
    #Epsilon - Mean Equation
    eps0m<-matAll$epsilon[i]
    #Parms Mean
    parmsM<-as.numeric(matAll[i,3:ncol(matAll)])

    #Create the data.frame
    treina.df<-as.matrix(cbind(wtrain, Xtrain))
    valida.df<-as.matrix(cbind(wvalid, Xvalid))

    #Estimate the SVM
    svmTreina<- CSVRL1(ytrain, treina.df, C0m, eps0m, kernel, parmsM)

    #Predict using validation dataset
    yPred <- PredictedCSVRL1(svmTreina, treina.df, valida.df)

    #Calculate the MSE
    MSE <- sum((yPred-yvalid)^2)

    #Return the results
    res<-data.frame(matAll[i,],"MSE"=MSE)
    res
  }

  #Stop clusters
  stopCluster(cl)

  #Get the minum MSE
  parmsEff <- svrCausality[which(svrCausality$MSE==min(svrCausality$MSE))[1],1:(ncol(svrCausality)-1)]

  #Create the complete data.frame
  df0Treina <- as.matrix(cbind(rep(0,nrow(Xtrain)), Xtrain))
  df0Valida <- as.matrix(cbind(rep(0,nrow(Xvalid)), Xvalid))
  df0 <- rbind(df0Treina,df0Valida)

  #Estimate the SVM
  yFinal0 <- c(ytrain, yvalid)
  svmFinal0 <- CSVRL1(yFinal0, df0, as.numeric(parmsEff[1]), as.numeric(parmsEff[2]), kernel, as.numeric(parmsEff[3:length(parmsEff)]))
  yPred0 <- PredictedCSVRL1(svmFinal0, df0, df0)

  #Create the complete data.frame
  df1Treina <- as.matrix(cbind(rep(1,nrow(Xtrain)), Xtrain))
  df1Valida <- as.matrix(cbind(rep(1,nrow(Xvalid)), Xvalid))
  df1 <- rbind(df1Treina,df1Valida)

  #Estimate the SVM
  yFinal1 <- c(ytrain, yvalid)
  svmFinal1 <- CSVRL1(yFinal1, df1, as.numeric(parmsEff[1]), as.numeric(parmsEff[2]), kernel, as.numeric(parmsEff[3:length(parmsEff)]))
  yPred1 <- PredictedCSVRL1(svmFinal1, df1, df1)

  #Final database
  XX <- rbind(Xtrain,Xvalid)
  WW <- c(wtrain,wvalid)
  final.df <- data.frame("Y"=yFinal1,"Y1"=yPred1,"Y0"=yPred0,"X"=XX,"W"=WW,"Tau"=yPred1-yPred0)

  #List of results
  list<-list("SVR"=svrCausality,"Final"=final.df,"SVM"=parmsEff,"MSE"=min(svrCausality$MSE)[1],"Kernel"=kernel)

  return(svrCausality)
}

