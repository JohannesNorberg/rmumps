// 
//


typedef struct spmatrix_ {
    int n;
    int nz;
    int *irn;
    int *jcn;
    double *a;
} spmatrix;




//////////


spmatrix makeRandomSpMatrix(int n, float fillRate);
