context("kernelmatrixcomputation")
library(mlRFinance)

test_that("KernelMatrixCamputation test", {
        set.seed(2929)
        p<-5
        n<-5
        dados<-matrix(rnorm(p*n),ncol=p,nrow=n)
        Kmat<-KernelMatrixComputation(dados,"Gaussian",c(2.5))
        expected <- matrix(
                      c(1.0000000, 0.4455850, 0.4522078, 0.7084599, 0.5825071,
                        0.4455850, 1.0000000, 0.1606014, 0.3566747, 0.7834800,
                        0.4522078, 0.1606014, 1.0000000, 0.6859084, 0.3400515,
                        0.7084599, 0.3566747, 0.6859084, 1.0000000, 0.6778726,
                        0.5825071, 0.7834800, 0.3400515, 0.6778726, 1.0000000),
                        5,5)
        expect_equal(Kmat, expected, tolerance=1e-07)
})

print("[SUCESS] All tests passed for 'kernelmatrixcomputation'!")