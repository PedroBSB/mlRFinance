#' Execute the Hansen, P. R. (2005) Superior Predictive Ability Test
#'
#' @param Dmat A matrix with the relative performance variables (NxK)
#' @param bVec A vector with the relative performance benchmark (Nx1)
#' @param typeFunc Type of function to be used:
#' 0 - max(0,x)
#' 1 - x*1_{x >= -sqrt{(w^2/n)2log log n}}
#' 2 - x
#' @param B total of boostraps samples.
#' @param geomMean mean of the geometric distribution used to generate the block lengths
#' @param alpha Type-I Error (level)
#' @param k Number of false rejections assumed
#' @param gamma The false discovery proportion (FDP) parameter
#' @return Statistic of the test and p-value
#' @examples
#' add(1, 1)
#' add(10, 1)
hansen.spa <- function(Dmat,bVec,typeFunc=1,B=1000,geomMean=20,bandwidth=0.5, alpha=0.05, k=1, gamma=0.1) {
  #Step 0: Computes the performance of model k relative to the benchmark at time t.
  d.mat <- apply(Dmat,2,function(x)x-bVec)

  #Size of the Time series
  n<-length(bVec)

  #Step 1: Compute d.bar
  d.bar  <- colMeans(d.mat)

  #Find a consistent estimator of V(n^(1/2)*overline(d)_{k})
  ##Computes the variance
  gamma0 <- sqrt(apply(d.mat,2,var))

  ## Compute the Auto-covariance
  auto.cov<-apply(d.mat,2, function(x) acf(x, "covariance", plot = F)$acf)
  auto.cov<-auto.cov[-1,]
  auto.cov<-rbind(gamma0,auto.cov)

  ##Compute the Kernel matrix Politis and Romano 1994
  i<-as.matrix(seq(1,nrow(d.mat)-1))
  kernel <- c(0,apply(i,1,function(x) ((n-x)/n)*((1-bandwidth)^x)+(x/n)*((1-bandwidth)^(n-x))))

  ##Compute a consistent estimator
  omega <-apply(auto.cov,2,function(x) x*kernel)
  omega <- auto.cov[1,]+colSums(omega)*2

  # Computes the T.spa Statistic
  t.SPA <- max(max(sqrt(n)*d.bar/omega),0)

  #Apply the Hansen function
  if(typeFunc==0){
    gFunc <- apply(as.matrix(d.bar),1,function(x) max(x,0))
  }else if(typeFunc==1){
    gFunc <- rep(NA,length(d.bar) )
    for(i in 1:length(d.bar)){
      gFunc[i] <- d.bar[i]*ifelse(d.bar[i]>=-sqrt(((omega[i]^2)/n)*log(log(n))),1,0)
    }
  }else{
    gFunc <- d.bar
  }
  Z <- t(apply(d.mat,1,function(x) x - gFunc))

  #Stationary Bootstrap
  Zboot <- boot::tsboot(Z,statistic=colMeans, R=B,l=geomMean,sim="geom")$t*sqrt(n)

  #Boostrap T.SPA
  T.SPABoot <- t(apply(Zboot,1,function(x) x/omega))

  #For each time series:
  T.SPA <- apply(as.matrix(apply(T.SPABoot,1,max)),1,function(x) max(x,0))

  #P-value
  p.value <- mean(T.SPA>t.SPA)

  #False Discovery Proportion (FDP)
  fdp <- FDPControl(t.SPA, t(T.SPABoot), gamma, alpha)

  #FamilyWise Error Rate (FWER)
  fwer <- FWERkControl(t.SPA, t(T.SPABoot), k, alpha)

  #Return function
  return(list("Hansen's SPA statistic"=t.SPA,
              "P-value"=p.value,
              "FDP"=fdp,
              "FWERk"=fwer))
}


#' Execute the White, H. (2000) Reality Check Test
#'
#' @param Dmat A matrix with the relative performance variables (NxK)
#' @param bVec A vector with the relative performance benchmark (Nx1)
#' @param typeFunc Type of function to be used:
#' 0 - max(0,x)
#' 1 - x*1_{x >= -sqrt{(w^2/n)2log log n}}
#' 2 - x
#' @param B total of boostraps samples.
#' @param geomMean mean of the geometric distribution used to generate the block lengths
#' @param alpha Type-I Error (level)
#' @param k Number of false rejections assumed
#' @param gamma The false discovery proportion (FDP) parameter
#' @return Statistic of the test and p-value
#' @examples
#' add(1, 1)
#' add(10, 1)
white.spa <- function(Dmat,bVec,typeFunc=1,B=1000,geomMean=20,bandwidth=0.5, alpha=0.05, k=1, gamma=0.1) {
  #Step 0: Computes the performance of model k relative to the benchmark at time t.
  d.mat <- apply(Dmat,2,function(x)x-bVec)

  #Size of the Time series
  n<-length(bVec)

  #Step 1: Compute d.bar
  d.bar  <- colMeans(d.mat)

  ##White works with a non standardized statistic
  omega <-rep(1,ncol(Dmat))

  # Computes the T.spa Statistic
  t.SPA <- max(max(sqrt(n)*d.bar/omega),0)

  #Apply the Hansen function
  if(typeFunc==0){
    gFunc <- apply(as.matrix(d.bar),1,function(x) max(x,0))
  }else if(typeFunc==1){
    gFunc <- rep(NA,length(d.bar) )
    for(i in 1:length(d.bar)){
      gFunc[i] <- d.bar[i]*ifelse(d.bar[i]>=-sqrt(((omega[i]^2)/n)*log(log(n))),1,0)
    }
  }else{
    gFunc <- d.bar
  }
  Z <- t(apply(d.mat,1,function(x) x - gFunc))

  #Stationary Bootstrap
  Zboot <- boot::tsboot(Z,statistic=colMeans, R=B,l=geomMean,sim="geom")$t*sqrt(n)

  #Boostrap T.SPA
  T.SPABoot <- t(apply(Zboot,1,function(x) x/omega))

  #For each time series:
  T.SPA <- apply(as.matrix(apply(T.SPABoot,1,max)),1,function(x) max(x,0))

  #P-value
  p.value <- mean(T.SPA>t.SPA)

  #False Discovery Proportion (FDP)
  fdp <- FDPControl(t.SPA, t(T.SPABoot), gamma, alpha)

  #FamilyWise Error Rate (FWER)
  fwer <- FWERkControl(t.SPA, t(T.SPABoot), k, alpha)

  #Return function
  return(list("White's RC statistic"=t.SPA,
              "P-value"=p.value,
              "FDP"=fdp,
              "FWERk"=fwer))
}



#' Execute the Classical Diebold and Mariano - EPA
#'
#' @param e.model loss function based on the comparasion
#' between model forecast and observed value (Nx1)
#' @param e.bench loss function based on the comparasion
#' between benchmark forecast and observed value (Nx1)
#' @param M lag considered
#' @return Statistic of the test and p-value
#' @examples
#' add(1, 1)
#' add(10, 1)
DM.epa <- function(e.model,e.bench, M){
  #If M is NA
  if(is.na(M)) M=ceiling(length(e.model)^(1/3))
  #Size of the serie
  T <- length(e.model)
  #Calculate the difference
  d.diff <- e.model-e.bench
  #Calculate the numerator
  d.bar <- mean(d.diff)
  #Compute the auto-covariance
  auto.cov <- acf(d.diff, "covariance", plot = F)$acf
  auto.cov <- auto.cov[2:M]
  #Calculate an estimate to spectral density
  var <- var(d.diff) +2*sum(auto.cov)
  #Diebold-Mariano Statistic
  DM <- d.bar/sqrt(var/T)
  #Calculate p-value
  pvalue <- 2*pnorm(-abs(DM))
  #Return function
  return(list("DM Statistic"=DM, "P-value"=pvalue))
}

#' Harvey, Leybourne, and Newbold (1997) Modification
#' to the Classical Diebold and Mariano - EPA
#'
#' @param e.model loss function based on the comparasion
#' between model forecast and observed value (Nx1)
#' @param e.bench loss function based on the comparasion
#' between benchmark forecast and observed value (Nx1)
#' @param M lag considered
#' @return Statistic of the test and p-value
#' @examples
#' add(1, 1)
#' add(10, 1)
#
DM.epa.corrected <- function(e.model,e.bench, M){
  if(is.na(M)) M=ceiling(length(e.model)^(1/3))
  #Size of the serie
  T <- length(e.model)
  #Calculate the difference
  d.diff <- e.model-e.bench
  #Calculate the numerator
  d.bar <- mean(d.diff)
  #Compute the auto-covariance
  auto.cov <- acf(d.diff, "covariance", plot = F)$acf
  auto.cov <- auto.cov[2:M]
  #Calculate an estimate to spectral density
  var <- var(d.diff) +2*sum(auto.cov)
  #Diebold-Mariano Statistic
  DM <- d.bar/sqrt(var/T)
  #Correct as Harvey, Leybourne, and Newbold (1997)
  DM <- DM*sqrt((T+1-2*M+M(M-1))/T)
  #Calculate p-value
  pvalue <-  2*pt(-abs(DM),df=T-1)
  #Return function
  return(list("Corrected DM Statistic"=DM, "P-value"=pvalue))
}
