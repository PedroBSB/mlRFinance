#Rodriguez-Lujan, Irene, et al. "Quadratic programming feature selection."
#Journal of Machine Learning Research 11.Apr (2010): 1491-1516.
QPFS<-function(y,X,alpha){
  #Compute the Similarity Matrix
  Q <- cor(X)
  Q <- abs(Q)
  #Compute the Relevance (binary problem)
  prob <- mean(y==+1)
  f <- diag(Q)*prob
  #Compute the weights
  weights <- QPFS(Q,f,alpha)
  #Return the results
  return(weights)
}

