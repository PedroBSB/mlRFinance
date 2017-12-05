#Page 160
A<-matrix(c(1,2,5,6),nrow=4,ncol=1)
d<-c(-1,-1,+1,-1)

svm2<- nuSVM(d, A, 0.5, "Polynomial", c(2,1))
svm2


