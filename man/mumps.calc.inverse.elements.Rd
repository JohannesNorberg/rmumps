\name{mumps.calc.inverse.elements}
\alias{mumps.calc.inverse.elements}
\title{Calculate elements of an inverse matrix.}
\description{
    Calculate chosen elements of an inverse of a sparse matrix.} 
\usage{
    mumps.calc.inverse.elements(mat,mask,sym = 0,np = detectCores())
}
\arguments{
    \item{mat}{Sparse matrix of class \code{\link[Matrix]{dgCMatrix-class}} or  \code{\link[Matrix]{dtCMatrix-class}}.}
    \item{mask}{Matrix of the same size as \code{mat}. Non-zero elements of \code{mask} are used as a mask for elements of the inverse matrix to be calculated.}
    \item{sym}{Symmetricity of matrix \code{mat}. \code{0} for unsymmetric, \code{1} for symmetric positive definite, or \code{2} for a general symmetric matrix.}
    \item{np}{Number of cores/threads used by MUMPS. }
}
\value{\code{mumps.calc.inverse.elements} returns a matrix of class 
\code{\link[Matrix]{dgCMatrix-class}} containing the elements of the inverse of \code{mat}.} 
\note{Matrix \code{mask} must be coercible to class \code{\link[Matrix]{dgCMatrix-class}}. }

\keyword{algebra,array}
