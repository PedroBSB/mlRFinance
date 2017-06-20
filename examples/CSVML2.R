#Page 160
A<-matrix(c(1,2,5,6),nrow=4,ncol=1)
d<-c(-1,-1,+1,-1)

svm2<- CSVML2(d, A, 50, "Polynomial", c(2,1),FALSE)
svm2





#Quadprog
K<-as.matrix(Matrix::nearPD(svm2$K)$mat)
Q<-as.matrix(Matrix::nearPD(svm2$Q)$mat)
matrixcalc::is.positive.definite(K)
dvec<-rep(-1,length(d))
I<-diag(1,nrow=length(d),ncol=length(d))
Amat<-rbind(I,-I)
C<-50
bvec <- c(rep(0,length(d)),rep(-C,length(d)))
quadprog::solve.QP(Q,-dvec,t(Amat),bvec,meq=0)

#Ipop 0.00000 14.90472 25.79428 16.26741
c<-rep(-1,length(d))
H<-Q
A<-I
b<-rep(0,length(d))
r<-rep(C,length(d))
l<-rep(0,length(d))
u<-rep(1e+10,length(d))
kernlab::ipop(c,H,A,b,l,u,r)
det(solve(H))



















library(kernlab)
X<-as.matrix(spam[,1:57])
y<-ifelse(spam$type=="nonspam",+1,-1)
svm2<- CSVML2(y, X, 50, "Polynomial", c(2,1))
svm2$SupportVectors<-round(svm2$SupportVectors,5)



#Quadprog
K<-apply(A,1,function(x) apply(A,1,function(y){
  (x*y+1)^2
}))
Q<-K*(d%*%t(d))

temp<-eigen(K)
temp$values[which(temp$values<0)]<-0

K.new<-temp$vectors%*%diag(temp$values)%*%t(temp$vector)


matrixcalc::is.positive.definite(as.matrix(K.new))

I<-diag(1,nrow=length(d),ncol=length(d))
Amat<-rbind(I,-I)

#Equality
Amat<-rbind(rep(1,length(d)),Amat)
C<-50
bvec <- c(1,rep(0,length(d)),rep(C,length(d)))
quadprog::solve.QP(Q,-dvec,t(Amat),bvec,meq=1)




################################################





svm1<- CSVML2(d, A, 1, "Gaussian", c(0.5))
svm1

svm2<- CSVML2(d, A, 1, "Polynomial", c(2,1))
svm2

#Quadprog
C<-1
sigma<-0.5
K<-apply(A,1,function(x) apply(A,1,function(y){
  exp(-sum((x-y)^2)/(2*sigma^2))
}))
K

D<-d%*%t(d)
CI<-diag(1,length(d))
Q<-(K+(CI/C))*D
dvec<-rep(1,length(d))

I<-diag(1,nrow=length(d),ncol=length(d))
Amat<-I
bvec <- c(rep(0,length(d)))

quadprog::solve.QP(Q,-dvec,Amat,-bvec,meq=0)



#Kernlab
library(kernlab)
svm<-ksvm(x=A,y=d,kernel ="rbfdot",kpar=list(sigma=1), type="C-svc", C=1, scaled = FALSE, cross=0)
svm@alpha[[1]]
svm1


