# Package: rmumps
# Simple (at the time being) R interface to MUMPS
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


read.dense.RHS <- function(filename) {

    fid <- file(filename,"rb")

    d1 <- readBin(fid,integer(),n=1,size=4)
    d2 <- readBin(fid,integer(),n=1,size=4)

    res <- readBin(fid,numeric(),n=d1*d2,size=8)

    dim(res) <- c(d1,d2)

    return(res)
}


mumps.solve <- function(mat,rhs,np = 4) {

    if (class(mat) != "dgCMatrix")
        stop("Matrix mat must be of class 'dgCMatrix'")

    # Write matrix and rhs to file
    save.sp.matrix(mat,"mumps_mat.bin")
    save.dense.RHS(rhs,"mumps_rhs.bin")

    # Run MUMPS

    # Construct command 
    # In this version sym=2 (general symmetric)
    cmd.path <- paste0(system.file(package="rmumps"),"/bin/mumpsdrv")
    command <- paste("mpirun","-np",np,cmd.path,"mumps_mat.bin","mumps_rhs.bin",0,sep=' ')
    # Execute command
    #cat(command,'\n')
    system(command)

    # Read solution from file
    sol <- read.dense.RHS("mumps_rhs.bin")
    # Construct matrix/vector
    #res <- matrix(sol$data,sol$rows,sol$cols)

    # Delete files

    # Return solution
    return(sol)
}
