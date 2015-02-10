# Package: rmumps
# Simple (for the time being) R interface to MUMPS
# Uses binary files and a executable MUMPS driver 


# Saves sparse matrix to binary file
save.sp.matrix <- function(mat,filename) {

    fid <- file(filename,"wb")

    # NB: This version assumes that matrix is square and symmetric
    n <- mat@Dim[1]
    nz <- length(mat@x)

    # save n
    writeBin(as.integer(n),fid,size = 4)
    # save nz
    writeBin(as.integer(nz),fid,size = 4)
    # save i
    ii <- mat@i + 1
    writeBin(as.integer(ii),fid,size=4)
    # save j
    jj <- rep(seq_along(diff(mat@p)),diff(mat@p))
    writeBin(as.integer(jj),fid,size=4)
    # save x
    writeBin(as.double(mat@x),fid,size=8)

    close(fid)
}


# Saves dense matrix to file
save.dense.RHS <- function(data,filename) {
    data <- as.matrix(data)
    dd <- dim(data)
    ll <- length(data)

    fid <- file(filename,"wb")

    writeBin(as.integer(dd[1]),fid,size=4)
    writeBin(as.integer(dd[2]),fid,size=4)

    writeBin(as.double(data[1:ll]),fid,size=8)

    close(fid)
}


# Reads dense matrix from file
read.dense.RHS <- function(filename) {

    fid <- file(filename,"rb")

    d1 <- readBin(fid,integer(),n=1,size=4)
    d2 <- readBin(fid,integer(),n=1,size=4)

    res <- readBin(fid,numeric(),n=d1*d2,size=8)

    dim(res) <- c(d1,d2)
    
    close(fid)

    return(res)
}


#######
# Main solve function
#
# Arguments:
#   mat     sparse matrix of class dgCMatrix or dtCMatrix
#   rhs     right hans side of the equation. Matrix or vector. Will be
#           saved as dense matrix.
#   np      Number of cores to be used in the calculation
#   sym     0 for non-symmetric, 1 for positive definite, 2 for
#           general symmetric
mumps.solve <- function(mat,rhs,np = 4,sym = 0) {

    # Check matrix class
    if (class(mat) != "dgCMatrix" && class(mat) !="dtCMatrix")
        stop("Matrix mat must be of class 'dgCMatrix' or 'dtCMatrix'")

    # If solving symmetric problem and matrix is given in dgCMatrix,
    # transform it into dtCMatrix format
    if ( sym > 0 && class(mat) == "dgCMatrix" ) {
        mat <- tril(mat)
    }


    # Write matrix and rhs to file
    save.sp.matrix(mat,"mumps_mat.bin")
    save.dense.RHS(rhs,"mumps_rhs.bin")

    # Run MUMPS

    # Construct command 
    mode <- 0
    cmd.path <- paste0(system.file(package="rmumps"),"/bin/mumpsdrv")
    command <- paste("mpirun","-np",np,cmd.path,mode,"mumps_mat.bin","mumps_rhs.bin",sym,sep=' ')
    # Execute command
    #cat(command,'\n')
    system(command)

    # Read solution from file
    sol <- read.dense.RHS("mumps_rhs.bin")
    # Construct matrix/vector
    #res <- matrix(sol$data,sol$rows,sol$cols)

    # Delete files
    file.remove("mumps_mat.bin")
    file.remove("mumps_rhs.bin")

    # Return solution
    return(sol)
}

mumps.calc.inverse.diagonal <- function(mat,np = 4, sym = 0) {
    # Check matrix class
    if (class(mat) != "dgCMatrix" && class(mat) !="dtCMatrix")
        stop("Matrix mat must be of class 'dgCMatrix' or 'dtCMatrix'")

    # If solving symmetric problem and matrix is given in dgCMatrix,
    # transform it into dtCMatrix format
    if ( sym > 0 && class(mat) == "dgCMatrix" ) {
        mat <- tril(mat)
    }


    # Write matrix and rhs to file
    save.sp.matrix(mat,"mumps_mat.bin")

    # Construct command 
    mode <- 1
    cmd.path <- paste0(system.file(package="rmumps"),"/bin/mumpsdrv")
    command <- paste("mpirun","-np",np,cmd.path,mode,"mumps_mat.bin",sym,sep=' ')
    # Execute command
    #cat(command,'\n')
    system(command)

    # Read diagonal values
    fid <- file("mumps_diag.bin","rb")
    
    n <- readBin(fid,integer(),1,size=4)

    res <- readBin(fid,numeric(),n,size=8)

    close(fid)

    file.remove("mumps_mat.bin")
    file.remove("mumps_diag.bin")

    return(res)


}
