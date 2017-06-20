#Aplication
#Page 160
A<-matrix(c(1,2,5,6),nrow=4,ncol=1)
d<-c(-1,-1,+1,-1)

K<-apply(A,1,function(x) apply(A,1,function(y){
  (x*y+1)^2
}))

#Q matrix
Q1<-K*(d%*%t(d))
D<-diag(d)
Q2<-D%*%K%*%D
all(Q1==Q2)


### Quad prog
library(quadprog)
Dmat<-Q1
dvec<-rep(1,nrow(A))
I<-diag(1,nrow=nrow(A),ncol=nrow(A))
Amat<-rbind(I,-I)
C<-50
bvec<-c(rep(0,nrow(A)),rep(-C,nrow(A)))

solve.QP(Dmat,dvec,t(Amat),bvec,meq=0)

### Ipop
library(kernlab)
c<-rep(-1,nrow(A))
H<-Q1
b<-rep(0,nrow(A))
r<-rep(C,nrow(A))
l<-rep(0,nrow(A))
u<-rep(1e+7,nrow(A))
A<-I
ipop(c, H, A, b, l, u, r)

help(solve.QP)

#Test
D<-diag(y)

Q1<-K*(y%*%t(y))



svm2<- CSVML1(y, X, 5, "Polynomial", c(2,1))
