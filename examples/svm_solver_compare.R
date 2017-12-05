#  This gist solves the hard-margin SVM problem in three ways: using quadprog, using kernlab's ipop, and by 
#  the e1071 wrapper around libsvm.
#
#  author: R. Walker  (r_walker@zoho.com)
#  LICENSE: MIT
library("quadprog")
library("kernlab")
library("e1071")

# Use Fisher iris data and binarize one of the species
# Choose "setosa" for a linearly separable example
# Choose "versicolor" or "virginica" for a non-separable example
# 
# For separable problem, all three methods should be very close.  For non-separable problem,
# ipop and quadprog should be close (after some parameter tuning) but SVM result will still be a bit
# different.
#
data(iris)
train <- iris
train$y <-ifelse(train[,5]=="setosa", 1, -1)

# order the training data labeling to avoid oddities with 
# training labels in libsvm
# see http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html#f430
train <- train[order(train$y, decreasing=TRUE),]

# set the problem data and parameters
X <- as.matrix(train[,c("Petal.Length", "Petal.Width")])
y <- as.matrix(train$y)
n <- dim(X)[1]

####################################################################################################################
# solve QP with quadprog and the perturbance hack
# From the documentation:
# This routine implements the dual method of Goldfarb and Idnani (1982, 1983) for solving quadratic programming
# problems of the form min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0.
####################################################################################################################
# quadprog solver requires that the D matrix be symmetric positive definite.  But the SVM problem is almost always
# only non-negative definite!  As a hack, we can perturb D by a small diagonal matrix and obtain positive definite
# matrix.  Choose eps a relatively small value for the diagonal perturbance.
eps <- 5e-4

# build the system matrices
Q <- sapply(1:n, function(i) y[i]*t(X)[,i])
D <- t(Q)%*%Q
d <- matrix(1, nrow=n)
b0 <- rbind( matrix(0, nrow=1, ncol=1) , matrix(0, nrow=n, ncol=1) )
A <- t(rbind(matrix(y, nrow=1, ncol=n), diag(nrow=n)))

# call the QP solver:
sol <- solve.QP(D +eps*diag(n), d, A, b0, meq=1, factorized=FALSE)
qpsol <- matrix(sol$solution, nrow=n)

####################################################################################################################
# solve the QP with kernlab solver
# ipop solves the quadratic programming problem :
# \min(c'*x + 1/2 * x' * H * x)
# subject to: 
# b <= A * x <= b + r
# l <= x <= u
#
# Unlike quadprog, ipop can handle a non-negative definite system matrix; however if the problem is 
# not separable then (as formulated here) the ipop solver will become unstable.  In this case, we can add
# the same perturbance term from above.
####################################################################################################################
uu <- 1e10     # Huge so that solution is unconstrained from above
eps <- 1e-11    # perturbance size (can set to 0 for separable problem)
b <- 0
r <- 0
A2 <- t(y)
l <- matrix(0, nrow=n, ncol=1)
u <- matrix(uu, nrow=n, ncol=1)
sol2 <- ipop(-d, t(Q)%*%Q+eps*diag(n), A2, b, l, u, r, verb=TRUE, sigf=5, margin=1e-8) # toggle margin and sigf if fails to converge.
ipopsol <- primal(sol2)
####################################################################################################################
# svm solver provided by e1071
# This solver is implemented to solve a soft-margin problem.  If we choose C very large, 
# then the margins will be very small and we'll approach the hard-margin classification problem that we solved above
# with ipop and quadprog
####################################################################################################################
C <- 1e5     # Huge value forces hard margin problem
sv <- svm(y~Petal.Length+Petal.Width, data=train, kernel="linear", scale=FALSE, type="C-classification", cost=C)

# get the slope and intercept terms from e1071 svm:
W <- rowSums(sapply(1:length(sv$coefs), function(i) sv$coefs[i]*sv$SV[i,]))
svmline = c(sv$rho/W[2], -W[1]/W[2])

# extract and plot the decision boundary from the results of each solver.
# build the support vectors, slopes, and intercepts
findLine <- function(a, y, X){
  nonzero <-  abs(a) > 1e-5
  W <- rowSums(sapply(which(nonzero), function(i) a[i]*y[i]*X[i,]))
  b <- mean(sapply(which(nonzero), function(i) X[i,]%*%W- y[i]))
  slope <- -W[1]/W[2]
  intercept <- b/W[2]
  return(c(intercept,slope))
}
qpline <- findLine( qpsol, y, X)
ipopline <- findLine(ipopsol, y, X)

# plot the results
library(ggplot2)
plt <- ggplot(train, aes(x=Petal.Length, y=Petal.Width)) + 
  ggtitle("Solving the SVM QP") +
  geom_point(aes(fill=factor(y)), size=3, pch=21) +
  geom_abline(intercept=qpline[1], slope=qpline[2], size=1, aes(color="quadprog"), show_guide=TRUE) +
  geom_abline(intercept=ipopline[1], slope=ipopline[2], size=1, aes(color="ipop")) +
  geom_abline(intercept=svmline[1], slope=svmline[2], size=1, aes(color="svm"))+
  scale_fill_manual(values=c("red","blue"), guide='none')+
  scale_color_manual(values=c("green", "yellow", "black"))+
  theme(legend.position="bottom", legend.title=element_blank())
print(plt)

# print the results
print(sprintf("quadprog number of nonzeros: %d", sum(abs(qpsol)>1e-7)))
print(sprintf("ipop number of nonzeros: %d", sum(abs(ipopsol)>1e-7)))
print(sprintf("svm number of nonzeros: %d", length(sv$coefs)))
print(sprintf("Difference between quadprog and ipop solutions: %f", sqrt(sum(qpsol - ipopsol)^2)))
print(sprintf("Quadprog: Intercept: %f    Slope: %f", qpline[1], qpline[2]))
print(sprintf("ipop: Intercept: %f    Slope: %f", ipopline[1], ipopline[2]))
print(sprintf("svm: Intercept: %f    Slope: %f", svmline[1], svmline[2]))