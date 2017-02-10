KPCA<-function(X,kernel,parms){
  #Get the Kernel Centered Matrix
  Kmod <- KPCAMatrix(X, kernel, parms)
  #Compute the PCA
  pca <- prcomp(Kmod,
                   center = FALSE,
                   scale = FALSE)
  #Return the results
  return(pca)
}

KPCApredict<-function(kpca,X,kernel,parms){
  #Get the Kernel Centered Matrix
  Kmod <- KPCAMatrix(X, kernel, parms)
  #Projection new data
  return(predict(kpca, newdata=Kmod))
}




