#Rodriguez-Lujan, Irene, et al. "Quadratic programming feature selection."
#Journal of Machine Learning Research 11.Apr (2010): 1491-1516.
#Only for binary problems
QPFS.binary<-function(y,X,alpha,type){
  ##### Compute the Similarity Matrix #####
  if(type=="Pearson"){
    Q <- cor(X)
    Q <- abs(Q)
  }
  ##### Compute the Relevance (binary problem) #####
  prob <- mean(y==+1)
  #First case (Case Y=+1)
  dummy <- rep(0,length(y))
  dummy[y==+1] <- +1
  #Estimate the correlation
  X.temp <- cbind(dummy,X)
  if(type=="Pearson"){
    Q.temp <- cor(X.temp)
    Q.temp <- abs(Q.temp)
  }
  f <- Q.temp[1,2:ncol(Q.temp)]*prob
  #Second case (Case Y=-1)
  prob <- mean(y==-1)
  #Dummy vector
  dummy <- rep(0,length(y))
  dummy[y==-1] <- +1
  #Estimate the correlation
  X.temp <- cbind(dummy,X)
  if(type=="Pearson"){
    Q.temp <- cor(X.temp)
    Q.temp <- abs(Q.temp)
  }
  f <- f + Q.temp[1,2:ncol(Q.temp)]*prob


  ##### Compute the weights #####
  weights <- QPFS(Q,f,alpha)
  #Return the results
  return(weights)
}

# X<-cbind(rnorm(1000),rnorm(1000),rnorm(1000),rnorm(1000))
# y<-rbinom(1000,1,0.2)
# y[y==0]<- -1
# alpha<-NA
# QPFS.binary(y,X,alpha,"Pearson")
