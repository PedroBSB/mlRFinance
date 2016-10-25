library(mlRFinance)
A<-matrix(0,ncol=3,nrow = 3)
A[1,1]<-A[2,2]<-A[3,3]<-1
A[1,2]<-A[2,1]<-3
A[1,3]<-A[3,1]<-2
A[2,3]<-A[3,2]<-4


nearPDefinite(A)

