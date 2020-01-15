
#include "henon.h"


void getStability(SStab M){
    int flag, dim;
    gsl_matrix* auxMat;
    gsl_eigen_nonsymmv_workspace * w;
    gsl_complex lambda; double realpart;

    // Init workspace
    dim = M.n;
    w = gsl_eigen_nonsymmv_alloc (dim);
    gsl_eigen_nonsymmv_params(1 , w);
    auxMat = gsl_matrix_alloc (dim, dim);

    // Copy matrix into it'c copy (needed since the solver destroys the matrix)
    flag = gsl_matrix_memcpy(auxMat, M.M);

    // Calculate eigenvectors and eigenvalues
    flag = gsl_eigen_nonsymmv(auxMat, M.eval, M.evec, w);
    gsl_eigen_nonsymmv_free (w);

    // Sort eigenvalues
    for (int i = 0; i < dim; i++){
        lambda = gsl_vector_complex_get (M.eval, i);
        realpart = gsl_complex_abs(lambda);
        if (realpart < 1.0) M.stable[i] = 1;
        else if (realpart > 1.0) M.stable[i] = 2;
        else if (fabs(realpart - 1.0) < 1.e-10) M.stable[i] = 3;
        else M.stable[i] = -1;
    }
}
