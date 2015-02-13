

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include "mpi.h"
#include "dmumps_c.h"

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I) - 1]



// Centralized assembled matrix
typedef struct spmatrix_ {
    int n;
    int nz;
    int *irn;
    int *jcn;
    double *a;
} spmatrix;

// Centralized dense right-hand side matrix
typedef struct densematrix_ {
    int nrhs;
    int lrhs;
    double *data;
} densematrix;

// Sparse right-hand side matrix
typedef struct elemspmatrix_ {
    int nz;
    int n;
    int *i;
    int *p;
    double *data;
} elemspmatrix;



spmatrix readSpMatrixFromFile(char *filename) {

    FILE *fid = fopen(filename,"rb");
    spmatrix mat;

    // Read matrix size
    fread(&mat.n,sizeof(int),1,fid);
    // Read nz
    fread(&mat.nz,sizeof(int),1,fid);

    // Allocate matrix data vectors
    mat.irn = malloc(mat.nz * sizeof(int));
    mat.jcn = malloc(mat.nz * sizeof(int));
    mat.a = malloc(mat.nz * sizeof(double));

    // Read data from file
    fread(mat.irn,sizeof(int),mat.nz,fid);
    fread(mat.jcn,sizeof(int),mat.nz,fid);
    fread(mat.a,sizeof(double),mat.nz,fid);

    fclose(fid);

    return mat;

}


densematrix readDenseMatrixFromFile(char *filename){

    FILE *fid = fopen(filename,"rb");
    densematrix mat;

    // Read RHS dimensions
    fread(&mat.lrhs,sizeof(int),1,fid);
    fread(&mat.nrhs,sizeof(int),1,fid);

    // Allocate data vector
    mat.data = malloc( mat.nrhs * mat.lrhs * sizeof(double));

    // Read data
    fread(mat.data,sizeof(double),mat.nrhs * mat.lrhs,fid);

    fclose(fid);

    return mat;
}


elemspmatrix readElemSpMatrixFromFile(char *filename) {
    FILE *fid = fopen(filename,"rb");
    elemspmatrix mat;

    // Read
    fread(&mat.nz,sizeof(int),1,fid);
    fread(&mat.n,sizeof(int),1,fid);

    mat.i = malloc(mat.nz * sizeof(int));
    mat.p = malloc((mat.n+1) * sizeof(int));
    mat.data = malloc(mat.nz * sizeof(double));

    fread(mat.i,sizeof(int),mat.nz,fid);
    fread(mat.p,sizeof(int),mat.n+1,fid);
    fread(mat.data,sizeof(double),mat.nz,fid);

    fclose(fid);

    return mat;
}



void writeDenseMatrixToFile(densematrix data, char *filename) {
    
    FILE *fid = fopen(filename,"wb");

    fwrite(&data.lrhs,sizeof(int),1,fid);
    fwrite(&data.nrhs,sizeof(int),1,fid);
    fwrite(data.data,sizeof(double),data.nrhs*data.lrhs,fid);

    fclose(fid);
}


void writeElemSpMatrixToFile(elemspmatrix mat, char *filename) {

    FILE *fid = fopen(filename,"wb");

    fwrite(&mat.nz,sizeof(int),1,fid);
    fwrite(&mat.n,sizeof(int),1,fid);
    fwrite(mat.i,sizeof(int),mat.nz,fid);
    fwrite(mat.p,sizeof(int),mat.n+1,fid);
    fwrite(mat.data,sizeof(double),mat.nz,fid);

    fclose(fid);
}


void mumps_solve(char *filename_mat, char *filename_rhs, int sym) {

    spmatrix mat;
    densematrix rhs;

    mat = readSpMatrixFromFile(filename_mat);
    rhs = readDenseMatrixFromFile(filename_rhs);

    DMUMPS_STRUC_C id;

    int myid, ierr;
    
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    id.job = JOB_INIT;
    id.par = 1;
    id.sym = sym;
    id.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&id);

    if (myid == 0) {
        id.n = mat.n;
        id.nz = mat.nz;
        id.irn = mat.irn;
        id.jcn = mat.jcn;
        id.a = mat.a;
        id.nrhs = rhs.nrhs;
        id.lrhs = rhs.lrhs;
        id.rhs = rhs.data;
    }

    id.ICNTL(3) = -1;
    id.ICNTL(4) = 1;

    //If sym==0, MATIS hangs the solver (in OSX)!
    if (sym ==0) {
        id.ICNTL(7) = 0;
    }
    else {
        id.ICNTL(7) = 7;
    }

    id.job = 6;
    dmumps_c(&id);


    if (myid==0) {
        // Get the solution
        int i;
        for ( i = 0 ; i < rhs.lrhs*rhs.nrhs ; i++) {
            rhs.data[i] = id.rhs[i];
        }

        // Save solution to file
        writeDenseMatrixToFile(rhs,filename_rhs);
    }


    id.job = JOB_END;
    dmumps_c(&id);

}


void mumps_diagonal(char *filename_mat, int sym) {

    spmatrix mat;

    mat = readSpMatrixFromFile(filename_mat);

    DMUMPS_STRUC_C id;

    int myid, ierr;
    
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    id.job = JOB_INIT;
    id.par = 1;
    id.sym = sym;
    id.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&id);

    int *ivec = malloc(mat.n * sizeof(int));
    int *pvec = malloc((mat.n + 1) * sizeof(int));
    double *dummyvec = malloc(mat.n * sizeof(double));
    int i;

    if (myid == 0) {

        // Construct sparse rhs representing the diagonal
        for (i = 0; i < mat.n ; i++) {
            ivec[i] = i+1;
            pvec[i] = i+1;
        }
        pvec[mat.n] = mat.n+1;

        id.n = mat.n;
        id.nz = mat.nz;
        id.irn = mat.irn;
        id.jcn = mat.jcn;
        id.a = mat.a;
        id.nz_rhs = mat.n;
        id.nrhs = mat.n;
        id.rhs_sparse = dummyvec;
        id.irhs_sparse = ivec;
        id.irhs_ptr = pvec;
    }
    

    id.ICNTL(3) = -1;
    id.ICNTL(4) = 1; 

    //If sym==0, MATIS hangs the solver!
    if (sym ==0) {
        id.ICNTL(7) = 0;
    }
    else {
        id.ICNTL(7) = 7;
    }

    // Calculate inverse matrix elements
    id.ICNTL(30) = 1;

    id.job = 6;
    dmumps_c(&id);

    if (myid==0) {
        // Get the solution
        int j;
        for ( j = 0 ; j < mat.n ; j++) {
            dummyvec[j] = id.rhs_sparse[j];
        }

        // Save diagonal to file
        FILE *fid = fopen("mumps_diag.bin","wb");

        fwrite(&mat.n,sizeof(int),1,fid);
        fwrite(dummyvec,sizeof(double),mat.n,fid);
        fclose(fid);
    }

    id.job = JOB_END;
    dmumps_c(&id);
}


void mumps_elements(char *filename_mat, char *filename_mask, int sym) {

    spmatrix mat;
    elemspmatrix mask;
    
    // Load matrices
    mat = readSpMatrixFromFile(filename_mat);
    mask = readElemSpMatrixFromFile(filename_mask);

    DMUMPS_STRUC_C id;

    int i;
    int myid, ierr;
    
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    id.job = JOB_INIT;
    id.par = 1;
    id.sym = sym;
    id.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&id);


    if (myid == 0) {
        id.n = mat.n;
        id.nz = mat.nz;
        id.irn = mat.irn;
        id.jcn = mat.jcn;
        id.a = mat.a;
        id.nz_rhs = mask.nz;
        id.nrhs = mask.n;
        id.rhs_sparse = mask.data;
        id.irhs_sparse = mask.i;
        id.irhs_ptr = mask.p;
    }
    

    id.ICNTL(3) = -1;
    id.ICNTL(4) = 1;

    //If sym==0, MATIS hangs the solver!
    if (sym ==0) {
        id.ICNTL(7) = 0;
    }
    else {
        id.ICNTL(7) = 7;
    }

    // Calculate inverse matrix elements
    id.ICNTL(30) = 1;

    id.job = 6;
    dmumps_c(&id);

    if (myid==0) {
        // Get the solution
        int j;
        for ( j = 0 ; j < mat.n ; j++) {
            mask.data[j] = id.rhs_sparse[j];
        }

        // Save diagonal to file
        writeElemSpMatrixToFile(mask,filename_mask);
    }

    id.job = JOB_END;
    dmumps_c(&id);

}



// Main routine
// Modes:
//  0 normal solve
//  1 calculate inverse matrix diagonal
//  2 calculate inverse matrix elements
//
int main(int argc, char *argv[]) {

    // Process arguments
    // First argument is the mode
    int mode = atoi(argv[1]);

    int ierr;
    ierr = MPI_Init(&argc,&argv);

    // Process remaining arguments based on mode
    switch(mode) 
    {
        case 0 :
            mumps_solve(argv[2],argv[3],atoi(argv[4]));
            break;
        case 1 :
            mumps_diagonal(argv[2],atoi(argv[3]));
            break;
        case 2 :
            mumps_elements(argv[2],argv[3],atoi(argv[4]));
            break;
        default :
            printf("Mode not implemented\n");
            exit(1);
    }

    MPI_Finalize();
    return 0;
}
