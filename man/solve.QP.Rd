\name{solve.QP}
\alias{solve.QP}
\title{Solve a Quadratic Programming Problem}
\description{
  This routine implements the dual method of Goldfarb and Idnani (1982,
  1983) for solving quadratic programming problems of the form
  \eqn{\min(-d^T b + 1/2 b^T D b)}{min(-d^T b + 1/2 b^T D b)} with the
  constraints \eqn{A^T b >= b_0}.
}
\usage{
solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE)
}
\arguments{
  \item{Dmat}{
    matrix appearing in the quadratic function to be minimized.
  }
  \item{dvec}{
    vector appearing in the quadratic function to be minimized.
  }
  \item{Amat}{
    matrix defining the constraints under which we want to minimize the
    quadratic function.
  }
  \item{bvec}{
    vector holding the values of \eqn{b_0} (defaults to zero).
  }
  \item{meq}{
    the first \code{meq} constraints are treated as equality
    constraints, all further as inequality constraints (defaults to 0).
  }
  \item{factorized}{
    logical flag: if \code{TRUE}, then we are passing
    \eqn{R^{-1}}{R^(-1)} (where \eqn{D = R^T R}) instead of the matrix
    \eqn{D}  in the argument \code{Dmat}.
  }
}
\value{
  a list with the following components:
  \item{solution}{
    vector containing the solution of the quadratic programming problem.
  }
  \item{value}{
    scalar, the value of the quadratic function at the solution
  }
  \item{unconstrained.solution}{
    vector containing the unconstrained minimizer of the quadratic
    function.
  }
  \item{iterations}{
    vector of length 2, the first component contains the number of
    iterations the algorithm needed, the second indicates how often
    constraints became inactive after becoming active first.
  }
  \item{Lagrangian}{
    vector with the Lagragian at the solution. 
  }
  \item{iact}{
    vector with the indices of the active constraints at the solution. 
  }
}

\references{
  D. Goldfarb and A. Idnani (1982).
  Dual and Primal-Dual Methods for Solving Strictly Convex
  Quadratic Programs.
  In J. P. Hennart (ed.), Numerical Analysis, Springer-Verlag,
  Berlin, pages 226--239.

  D. Goldfarb and A. Idnani (1983).
  A numerically stable dual method for solving strictly convex quadratic 
  programs.
  \emph{Mathematical Programming}, \bold{27}, 1--33.
}  

\seealso{
  \code{\link{solve.QP.compact}}
}

\examples{
##
## Assume we want to minimize: -(0 5 0) \%*\% b + 1/2 b^T b
## under the constraints:      A^T b >= b0
## with b0 = (-8,2,0)^T
## and      (-4  2  0) 
##      A = (-3  1 -2)
##          ( 0  0  1)
## we can use solve.QP as follows:
##
Dmat       <- matrix(0,3,3)
diag(Dmat) <- 1
dvec       <- c(0,5,0)
Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)
bvec       <- c(-8,2,0)
solve.QP(Dmat,dvec,Amat,bvec=bvec)
}

\keyword{optimize}

