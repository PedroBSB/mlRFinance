ipoptr_qp <- function(Dmat, dvec, Amat, bvec, ub=100){
  # Solve the quadratic program
  #
  # min -d^T x + 1/2 x^T D x
  #   s.t. A%*%x>= b
  #
  # with ipoptr.
  
  n <- length(bvec)
  # Jacobian structural components
  j_vals <- unlist(Amat[Amat!=0])
  j_mask <- make.sparse(ifelse(Amat!=0, 1, 0))

  # Hessian structural components
  h_mask <- make.sparse(ifelse(upper.tri(Dmat, diag=TRUE) & Dmat!=0, 1, 0))
  h_vals <- do.call(c, sapply(1:length(h_mask), function(i) Dmat[i,h_mask[[i]]]))

  # build the ipoptr inputs
  eval_f <- function(x) return(-t(dvec)%*%x + 0.5*t(x)%*%Dmat%*%x)
  eval_grad_f <- function(x) return(-dvec + Dmat%*%x)
  eval_h <- function(x, obj_factor, hessian_lambda) return(obj_factor*h_vals)
  eval_h_structure <- h_mask
  eval_g <- function(x) return(t(Amat)%*%x)
  eval_jac_g <- function(x) return(j_vals)
  eval_jac_g_structure <- j_mask
  constraint_lb <- bvec
  constraint_ub <- rep( ub, n)

  # initialize with the global unconstrained minimum, as done in quadprog package
  # NOTE: This will only work if lb <= x0 <= ub.  If this is not the case, 
  # use x0 = lb can be used instead.
  x0 <- solve(Dmat, dvec)
    
  # call the solver
  sol <- ipoptr(x0 = x0, 
  eval_f = eval_f,
  eval_grad_f = eval_grad_f, 
  eval_g = eval_g, 
  eval_jac_g = eval_jac_g,
  eval_jac_g_structure = eval_jac_g_structure,
  constraint_lb = constraint_lb,
  constraint_ub = constraint_ub,
  eval_h = eval_h,
  eval_h_structure = eval_h_structure)
  return(sol)
}