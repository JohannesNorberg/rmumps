

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "spmatrix.h"







spmatrix readSpMatrixFromFile(char *filename) {

    FILE *fid = fopen(filename,"rb");
    spmatrix mat;
    int status;

    // Read matrix size
    fread(&mat.n,sizeof(int),1,fid);
    // Read nz
    fread(&mat.nz,sizeof(int),1,fid);

    // Allocate matrix data vectors
    mat.irn = malloc(mat.nz * sizeof(int));
    mat.jcn = malloc(mat.nz * sizeof(int));
    mat.a = malloc(mat.nz * sizeof(double));

    // Read data from file
    status = fread(mat.irn,sizeof(int),mat.nz,fid);
    printf("ss: %d\n",status);
    status = fread(mat.jcn,sizeof(int),mat.nz,fid);
    printf("ss: %d\n",status);
    status = fread(mat.a,sizeof(double),mat.nz,fid);
    printf("ss: %d\n",status);

    fclose(fid);

    return mat;

}




int main(int argc, char *argv[]) {

    char *filename;
    spmatrix mat;
    
    filename = argv[1];

    printf("Reading file %s\n",filename);

    mat = readSpMatrixFromFile(filename);

    printf("n: %d\n",mat.n);
    printf("nz: %d\n",mat.nz);
   
    int i;
    for (i = 0 ; i < 20 ; i++) {
        printf("%d: %d %d %lf\n",i,mat.irn[i],mat.jcn[i],mat.a[i]);
    }

    return 0;

}
