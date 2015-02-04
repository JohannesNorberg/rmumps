// Construct sparse matrices matrices for MUMS
// Used for testing purposes
//

#include<math.h>
#include<stdlib.h>
#include<time.h>
#include "spmatrix.h"




spmatrix makeRandomSpMatrix(int n, float fillRate) {
    // Init spmatrix
    spmatrix mat;
    mat.n = n;

    int nz = (int) (fillRate * n * n ), rr, ss;

    mat.irn = malloc(nz * sizeof(int));
    mat.jcn = malloc(nz * sizeof(int));
    mat.a = malloc(sizeof(double));
    mat.nz = nz;

    // Fill the matrix (with ones)
    // Full diagonal

    int i;
    for ( i = 0 ; i < n ; i++) {
        mat.irn[i] = i;
        mat.jcn[i] = i;
        mat.a[i] = 1.0f;
    }
    
    //time_t t;
    //srand((unsigned) time(&t));

    // Fill the rest of the matrix
    if (nz - n > 0) {
        for (i = n; i < nz ; i++) {
            rr = rand() % n;
            ss = rand() % n;
            mat.irn[i] = rr;
            mat.jcn[i] = ss;
            mat.a[i] = 1.0f;
    }

    return mat;
}



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
    fread(mat.irn,mat.nz * sizeof(int),mat.nz,fid);
    fread(mat.jcn,mat.nz * sizeof(int),mat.nz,fid);
    fread(mat.a,mat.nz * sizeof(double),mat.nz,fid);

    close(fid);

    return mat;

}
