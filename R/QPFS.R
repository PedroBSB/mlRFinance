#Rodriguez-Lujan, Irene, et al. "Quadratic programming feature selection."
#Journal of Machine Learning Research 11.Apr (2010): 1491-1516.
QPFS<-function(y,X,alpha,type){
  #Compute the Similarity Matrix
  if(type="Pearson"){
    Q <- cor(X)
    Q <- abs(Q)
  }
  #Compute the Relevance (binary problem)
  prob <- mean(y==+1)
  f <- diag(Q)*prob
  #Compute the weights
  weights <- QPFS(Q,f,alpha)
  #Return the results
  return(weights)
}

# X<-cbind(rnorm(1000),rnorm(1000),rnorm(1000),rnorm(1000))
# y<-rbinom(1000,1,0.2)
# y[y==0]<- -1
# alpha<-NA
