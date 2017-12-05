context("quadprog")
library(mlRFinance)

test_that("quaprogInequality test1", {
    sol<-c(0, 0, 0)
    Dmat       <- matrix(0,3,3)
    diag(Dmat) <- 1
    dvec       <- c(0,5,0)
    Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)
    bvec       <- c(-8,2,0)
    sol2<-quadprogInequality(sol, Dmat, dvec, t(Amat), -bvec)
    expect_equal(sol2, matrix(c(-0.4651163, -1.0697674, -2.1395349), 3, 1), tolerance=1e-07)
})

test_that("quadprogInequality test2", {
    Dmat <- matrix(0,2,2)
    diag(Dmat) <- 2
    dvec <-c(-8,-6)
    Amat <- matrix(c(-1, 0,0,-1,1, 1),3,2,byrow=TRUE)
    bvec <- c(0,0,5)
    sol<-c(0,0)
    sol1<-quadprogInequality(sol, Dmat, dvec, Amat, bvec)
    expect_equal(sol1, matrix(c(3, 2),2,1), tolerance=1e-07)
})

test_that("quadprogInequality test3",{
    sol<-c(0, 0)
    Dmat       <- matrix(0,2,2)
    Dmat[1,1]<-2
    Dmat[2,2]<-8
    dvec       <- c(-8,-16)
    Amat       <- matrix(c(1, 1,1, 0,-1, 0,0,-1),4,2,byrow=T)
    bvec       <- c(5,3,0,0)
    sol3<-quadprogInequality(sol, Dmat, dvec, Amat, bvec)
    expect_equal(sol3, matrix(c(3, 2),2,1), tolerance=1e-07)
})


print("[SUCESS] All tests passed for 'quadprog'!")
