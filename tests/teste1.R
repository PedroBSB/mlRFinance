#mlRFinance - Example 1
sol<-c(0, 0, 0)
Dmat       <- matrix(0,3,3)
diag(Dmat) <- 1
dvec       <- c(0,5,0)
Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)
bvec       <- c(-8,2,0)
sol2<-quadprogInequality(sol, Dmat, dvec, t(Amat), -bvec)
sol2

#mlRFinance
CE <- matrix(0, nrow(Dmat), 0)
ce0 <- vector()
rcppeigen_quadratic_solve(Dmat,-dvec, CE, -ce0, Amat, -bvec)
solveTest(Dmat, dvec, Amat, bvec, CE,  ce0 )
solveTest2(Dmat, dvec, Amat, bvec )
#Quadprog
quadprog::solve.QP(Dmat, dvec, Amat, bvec)


#mlRFinance - Example 2
Dmat <- matrix(0,2,2)
diag(Dmat) <- 2
dvec <-c(-8,-6)
Amat <- matrix(c(-1, 0,
                 0,-1,
                 1, 1),3,2,byrow=TRUE)
bvec <- c(0,0,5)
sol<-c(0,0)
sol1<-quadprogInequality(sol, Dmat, dvec, Amat, bvec)
sol1
#Quadprog
quadprog::solve.QP(Dmat,dvec,t(Amat),-bvec)
#mlRFinance
solveTest2(Dmat, dvec, t(Amat),-bvec )



#mlRFinance - Example 3
sol<-c(0, 0)
Dmat       <- matrix(0,2,2)
Dmat[1,1]<-2
Dmat[2,2]<-8
dvec       <- c(-8,-16)
Amat       <- matrix(c(1, 1,
                       1, 0,
                       -1, 0,
                       0,-1),4,2,byrow=T)
bvec       <- c(5,3,0,0)
sol3<-quadprogInequality(sol, Dmat, dvec, Amat, bvec)
sol3
quadprog::solve.QP(Dmat,dvec,t(Amat),-bvec)
solveTest2(Dmat, dvec, t(Amat),-bvec )
