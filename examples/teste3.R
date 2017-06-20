G = matrix(c(2.1, 1.5, 1.2,
             0.0, 2.2, 1.3,
             1.0, 0.0, 3.1), 3, 3)

g0 = c(6.0, 1.0, 1.0)

CE = matrix(c(1, 2, -1), 3, 1)

ce0 = c(-4);

CI = matrix(c( 1,  0,  0,
               0,  1,  0,
               0,  0,  1,
               -1, -1,  0), 3, 4)

ci0 = c(0, 0, 0, 10)

rcppeigen_quadratic_solve(G,g0,CE,ce0,CI,ci0)
