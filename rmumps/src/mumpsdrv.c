

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include "mpi.h"
#include "dmumps_c.h"

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I) - 1]




typedef struct spmatrix_ {
    int n;
    int nz;
    int *irn;
    int *jcn;
    double *a;
} spmatrix;

typedef struct denserhs_ {
    int nrhs;
    int lrhs;
    double *data;
} denserhs;





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


denserhs readDenseRHSFromFile(char *filename){

    FILE *fid = fopen(filename,"rb");
    denserhs mat;

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



void writeDenseRHStoFile(denserhs data, char *filename) {
    
    FILE *fid = fopen(filename,"wb");

    fwrite(&data.lrhs,sizeof(int),1,fid);
    fwrite(&data.nrhs,sizeof(int),1,fid);

    fwrite(data.data,sizeof(double),data.nrhs*data.lrhs,fid);

    fclose(fid);
}




int main(int argc, char *argv[]) {

    spmatrix mat;
    denserhs rhs;
    char *filename_mat, *filename_rhs;

    filename_mat = argv[1];
    filename_rhs = argv[2];
    int sym = atoi(argv[3]);

//    printf("Parameters:\n %s\n %s\n %d   \n\n\n",filename_mat,filename_rhs,sym);

    mat = readSpMatrixFromFile(filename_mat);
    rhs = readDenseRHSFromFile(filename_rhs);

    //printf("mat: %f %f %f %f\n",mat.a[0],mat.a[1],mat.a[2],mat.a[3]);
    //printf("rhs: %f %f %f %f\n\n",rhs.data[0],rhs.data[1],rhs.data[2],rhs.data[3]);
//    printf("mat: %d %d\n",mat.n,mat.nz);
//    printf("rhs: %d %d\n",rhs.lrhs,rhs.nrhs);


    DMUMPS_STRUC_C id;

    int myid, ierr;
    ierr = MPI_Init(&argc,&argv);
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

    //id.ICNTL(1) = -1;
    //id.ICNTL(2) = -1;
    //id.ICNTL(3) = -1;
    id.ICNTL(4) = 2; // HOX! Change this when ready!!


    id.ICNTL(7) = 0;
    //id.ICNTL(28) = 2;
    //id.ICNTL(29) = 2;

    id.job = 6;
    dmumps_c(&id);


    if (myid==0) {
        // Get the solution
        int i;
        for ( i = 0 ; i < rhs.lrhs*rhs.nrhs ; i++) {
            rhs.data[i] = id.rhs[i];
        }

        // Save solution to file
        writeDenseRHStoFile(rhs,filename_rhs);
    }


    id.job = JOB_END;
    dmumps_c(&id);

    MPI_Finalize();

    return 0;

}
