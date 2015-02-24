\name{mumps.calc.inverse.diagonal}
\alias{mumps.calc.inverse.diagonal}
\title{Calculate diagonal of an inverse matrix.}
\description{
    Calculate the diagonal elements of an inverse of a sparse matrix, i.e. \code{diag(solve(mat))}.
    }
\usage{
    mumps.calc.inverse.diagonal(mat,sym = 0,np = detectCores())
}
\arguments{
    \item{mat}{Sparse matrix of class \code{\link[Matrix]{dgCMatrix-class}} or  \code{\link[Matrix]{dtCMatrix-class}}.}
    \item{sym}{Symmetricity of matrix \code{mat}. \code{0} for unsymmetric, \code{1} for symmetric positive definite, or \code{2} for a general symmetric matrix.}
    \item{np}{Number of cores/threads used by MUMPS. }
}
\value{\code{mumps.calc.inverse.diagonal} returns a vector containing the diagonal values of the inverse of \code{mat}.} 


\keyword{algebra,array}
