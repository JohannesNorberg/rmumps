\name{mumps.solve}
\alias{mumps.solve}
\title{Solve sparse equation}
\description{
    Solve \code{x} from a sparse matrix equation \code{mat\%*\%x = rhs} using MUMPS.
    }
\usage{
    mumps.solve(mat,rhs,sym = 0,np = detectCores())
}
\arguments{
   \item{mat}{Sparse matrix of class \code{\link[Matrix]{dgCMatrix-class}} or  \code{\link[Matrix]{dtCMatrix-class}}.}
    \item{rhs}{Vector or matrix containing the right hand side of the equation. Both sparse and dense matrices are supported.}
    \item{sym}{Symmetricity of matrix \code{mat}. \code{0} for unsymmetric, \code{1} for symmetric positive definite, or \code{2} for a general symmetric matrix.}
    \item{np}{Number of cores/threads used by MUMPS. }
}
\value{\code{mumps.solve} returns the solution to the equation as a dense vector/matrix. Nota that the solution is always returned as a standard R matrix of class \code{\link{matrix}}.}


\keyword{algebra,array}
