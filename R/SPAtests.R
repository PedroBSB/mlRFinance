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
#' @return Statistic of the test and p-value
#' @examples
#' add(1, 1)
#' add(10, 1)
hansen.spa <- function(Dmat,bVec,typeFunc,B,geomMean,bandwidth) {
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
    gFunc <- apply(as.matrix(d.bar),1,function(x) x*ifelse(x>=-sqrt(((omega^2)/n)*log(log(n))),1,0))
  }else{
    gFunc <- d.bar
  }
  Z <- t(apply(d.mat,1,function(x) x - gFunc))

  #Stationary Bootstrap
  Zboot <- boot::tsboot(Z,statistic=colMeans, R=B,l=geomMean,sim="geom")$t*sqrt(n)

  #Boostrap T.SPA
  T.SPA <- t(apply(Zboot,1,function(x) x/omega))

  #For each time series:
  T.SPA <- apply(as.matrix(apply(T.SPA,1,max)),1,function(x) max(x,0))

  #P-value
  p.value <- mean(T.SPA>t.SPA)

  #Return function
  return(list("Hansen's SPA statistic"=t.SPA,
              "P-value"=p.value))
}



