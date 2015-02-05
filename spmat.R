# Create a symmetric positive definite sparse matrix and save it to binary file
library(Matrix)
make.sym.spmat <- function(n,fill,fname){
  
  #A <- sparseMatrix(i=1,j=1,x=1,dims = c(n,n))
  
  ii <- 1:n
  jj <- 1:n
  xx <- rep(1,n)
  
  nz <- round(fill*n**2)
  
  r1 <- ceiling(runif(nz-n,0,n))
  r2 <- ceiling(runif(nz-n,0,n))
  xx2 <- rep(1,nz-n)
  
  A <- sparseMatrix(i=c(ii,r1),j=c(jj,r2),x=c(xx,xx2),dims = c(n,n))
  
  
  
  cat('Matrix ready\n')
  
  B <- t(A)%*%A
  
  nz <- length(B@x)
  
  # Save matrix B in binary format
  fid <- file(fname,"wb")
  
  # save n
  writeBin(as.integer(n),fid,size = 4)
  # save nz
  writeBin(as.integer(nz),fid,size = 4)
  # save i
  ii <- B@i + 1
  writeBin(as.integer(ii),fid,size=4)
  # save j
  jj <- rep(seq_along(diff(B@p)),diff(B@p))
  writeBin(as.integer(jj),fid,size=4)
  # save x
  writeBin(as.double(B@x),fid,size=8)
  
  close(fid)
  
  return(B)
  
}