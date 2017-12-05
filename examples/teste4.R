#http://stackoverflow.com/questions/33088363/how-to-use-r-package-quadprog-to-solve-svm
library(kernlab) # for the spam data
# Load the input data to be used
data(spam)

# Use only a subset of the data (20%)
spam <- spam[sample(nrow(spam), round(0.2 * nrow(spam)), replace = FALSE), ]

# Retrieve the features and data
X <- spam[, 1:(ncol(spam) - 1)]
Y_f <- spam[, ncol(spam)]
Y <- 2 * (as.numeric(Y_f) - 1.5) # {-1, 1}

# Sample size
N <- nrow(X)
# Number of dimensions
n_d <- ncol(X)

# Value of the regularization parameter
C <- 1

d <- c(0, rep(0, n_d), rep(-C, N))

eps <- 1e-10
D <- diag(c(eps, rep(1, n_d), rep(eps, N)))

I_N <- diag(N) # N x N identity matrix

A_1 <- cbind(matrix(0, ncol = n_d + 1, nrow = N), I_N)
A_2 <- as.matrix(cbind(as.matrix(Y), X * as.matrix(Y)[, rep(1, n_d)], I_N))
rownames(A_1) <- NULL; rownames(A_2) <- NULL
colnames(A_1) <- NULL; colnames(A_2) <- NULL

A <- t(rbind(A_1, A_2))
b_0 <- c(rep(0, N), rep(1, N))

library(quadprog)
results <- solve.QP(D, d, A, b_0)

# Retrieve the results
b_optim <- results$solution



svm1<- CSVML1(d, A, C, "Gaussian", c(0.5))
