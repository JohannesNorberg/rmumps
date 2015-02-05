

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include "spmatrix.h"
#include "mpi.h"
#include "dmumps_c.h"

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I) - 1]


int main(int argc, char *argv[]) {

    printf("Starting...\n");
    spmatrix mat;
    char *filename;

    filename = argv[1];
    int ss = atoi(argv[2]);

    printf(" Reading data from file...");
    mat = readSpMatrixFromFile(filename);
    printf("DONE\n");

    printf(" Making rhs...");
    double *rhs = malloc(sizeof(double) * mat.n);

    int i;
    for (i = 0; i < mat.n ; i++) rhs[i] = 1.0f;
    printf("DONE\n");

    DMUMPS_STRUC_C id;

    printf(" Initializing MPI...");
    int myid, ierr;
    ierr = MPI_Init(&argc,&argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    printf("DONE\n");

    id.job = JOB_INIT;
    id.par = 1;
    id.sym = ss;
    id.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&id);

    if (myid == 0) {
        id.n = mat.n;
        id.nz = mat.nz;
        id.irn = mat.irn;
        id.jcn = mat.jcn;
        id.a = mat.a;
        id.rhs = rhs;
    }

    //id.ICNTL(1) = -1;
    //id.ICNTL(2) = -1;
    //id.ICNTL(3) = -1;
    id.ICNTL(4) = 2;

    id.job = 6;
    dmumps_c(&id);

    id.job = JOB_END;
    dmumps_c(&id);

    if(myid == 0) {
        printf("%f %f %f\n",id.rhs[0],id.rhs[5000],id.rhs[9999]);
    }

    return 0;

}
