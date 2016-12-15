#Page 160 - Small
library(mlRFinance)
A<-matrix(c(1,2,5,6),nrow=4,ncol=1)
d<-c(-1,-1,+1,-1)

svm2<- CSVML1(d, A, 50, "Polynomial", c(2,1), FALSE)
svm2

R2PredictedCSVML1(svm2,d, A, 2, FALSE)


K<-as.matrix(Matrix::nearPD(svm2$K)$mat)
Q<-as.matrix(Matrix::nearPD(svm2$Q)$mat)

#Quadprog
matrixcalc::is.positive.definite(K)
dvec<-rep(-1,length(d))
I<-diag(1,nrow=length(d),ncol=length(d))
Amat<-rbind(I,-I)
C<-50
bvec <- c(rep(0,length(d)),rep(-C,length(d)))
quadprog::solve.QP(Q,-dvec,t(Amat),bvec,meq=0)

#Ipop
c<-rep(-1,length(d))
H<-Q
A<-I
b<-rep(0,length(d))
r<-rep(C,length(d))
l<-rep(0,length(d))
u<-rep(1e+10,length(d))
kernlab::ipop(c,H,A,b,l,u,r)
det(solve(H))

#Aplication

library(kernlab)
data(spam)
X<-as.matrix(spam[,1:57])
y<-ifelse(spam$type=="nonspam",+1,-1)
svm2<-rep(0,length(y))

# Start the clock!
ptm <- proc.time()

svm2<- CSVML1(y, X, 50, "Polynomial", c(2,1), FALSE)

# Stop the clock
proc.time() - ptm


#Positive Definite
Q<-as.matrix(Matrix::nearPD(svm2$Q)$mat)

#Quadprog
dvec<-rep(-1,length(y))
I<-diag(1,nrow=length(y),ncol=length(y))
Amat<-rbind(I,-I)

#Equality
C<-50
bvec <- c(rep(0,length(y)),rep(-C,length(y)))
quadprog::solve.QP(Q,-dvec,t(Amat),bvec,meq=0)


#Ipop
help(ipop)
c<-rep(-1,length(y))
H<-Q
A<-I
b<-rep(0,length(y))
r<-rep(C,length(y))
l<-rep(0,length(y))
u<-rep(1e+10,length(y))
kernlab::ipop(c,H,A,b,l,u,r)


#########################################################








#svm1<- CSVML1(d, A, 1, "Gaussian", c(0.5), FALSE)
#svm1

#Quadprog
sigma<-0.5
K<-apply(A,1,function(x) apply(A,1,function(y){
  exp(-sum((x-y)^2)/(2*sigma^2))
}))
Q<-diag(d)%*%K%*%diag(d)
dvec<-rep(1,length(d))

I<-diag(1,nrow=length(d),ncol=length(d))
Amat<-rbind(I,-I)

#Equality
Amat<-rbind(rep(1,length(d)),Amat)
C<-10
bvec <- c(1,rep(0,length(d)),rep(C,length(d)))
quadprog::solve.QP(Q,dvec,t(Amat),bvec,meq=1)



#Kernlab
library(kernlab)
svm<-ksvm(x=A,y=d,kernel ="rbfdot",kpar=list(sigma=0.5), type="C-svc", C=1, scaled = FALSE, cross=0)
svm@alpha[[1]]

