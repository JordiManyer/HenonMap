
#include "henon.h"


void getParametrization(Param parStable , Param parUnst , int dim , Point x0 , SStab M , double eps) {
    int iStable, iUnst;
    double eigStable, eigUnst;
    gsl_complex lambda; double realpart;
    gsl_vector* v0; double aux;

    // Get stable and unstable eigenvalues and their positions
    iStable = -1; iUnst = -1;
    eigStable = 100.0; eigUnst = 0.0;
    for (int i = 0; i < dim ; ++i) {
        lambda = gsl_vector_complex_get (M.eval, i);
        if (GSL_IMAG(lambda) > 1.e-10) printf("WARNING: Eigenvalue is not full real!!");
        realpart = GSL_REAL(lambda);

        if (fabs(realpart) < fabs(eigStable)) {iStable = i; eigStable = realpart;}
        if (fabs(realpart) > fabs(eigUnst)) {iUnst = i; eigUnst = realpart;}
    }
    printf("\ngetParametrization.c : \n");
    printf("iStable: %d , iUnstable: %d \n" , iStable , iUnst);
    printf("eigStable: %g , eigUnstable: %g \n" , eigStable , eigUnst);

    // UNSTABLE MANIFOLD:
    v0 = gsl_vector_calloc(dim);
    for (int i = 0; i < dim; ++i) {
        if (GSL_IMAG(gsl_matrix_complex_get(M.evec , iUnst , i)) > 1.e-10) printf("WARNING: Eigenvector is not full real!!");
        gsl_vector_set(v0 , i , GSL_REAL(gsl_matrix_complex_get(M.evec , iUnst , i)));
    }

    paramMethod1D( parUnst , dim , x0 , eigStable , v0 , eps);
    gsl_vector_free(v0);

    // STABLE MANIFOLD:
    v0 = gsl_vector_calloc(dim);
    for (int i = 0; i < dim; ++i) {
        if (GSL_IMAG(gsl_matrix_complex_get(M.evec , iStable , i)) > 1.e-10) printf("WARNING: Eigenvector is not full real!!");
        gsl_vector_set(v0 , i , GSL_REAL(gsl_matrix_complex_get(M.evec , iStable , i)));
    }

    paramMethod1D( parStable , dim , x0 , eigUnst , v0 , eps);
    gsl_vector_free(v0);
}



void paramMethod1D(Param par , int dim , Point p0 , double eval , gsl_vector* evec , double eps) {
    int flag ; int sign;
    gsl_matrix* A; gsl_vector* b; gsl_permutation* perm;
    gsl_vector* aux; double sum;

    // Initial values
    par.x[0] = p0.x; par.y[0] = p0.y;
    par.x[1] = gsl_vector_get(evec , 0); par.y[1] = gsl_vector_get(evec , 1);

    // Allocate system
    A = gsl_matrix_calloc(dim , dim);
    b = gsl_vector_calloc(dim);
    perm = gsl_permutation_calloc(dim);
    aux = gsl_vector_calloc(dim);

    // Loop to calculate the i-th coefficients
    for (int i = 2 ; i < par.n ; ++i) {

        // Calculate A
        gsl_matrix_set(A , 0 , 0 , 1.0 + eps*eps - 2.0*eps*eps*p0.x - pow(eval,i));
        gsl_matrix_set(A , 1 , 0 , eps - 2.0*eps*p0.x);
        gsl_matrix_set(A , 0 , 1 , eps);
        gsl_matrix_set(A , 1 , 1 , 1-pow(eval,i));

        // Calculate b
        sum = 0.0;
        for (int k = 1; k < i; ++k) {
            sum += par.x[i-k]*par.x[k];
        }
        gsl_vector_set(b , 0 , eps*eps*sum);
        gsl_vector_set(b , 1 , eps*sum);

        // Solve the system A x = b
        gsl_permutation_init(perm); // Init permutation
        flag = gsl_linalg_LU_decomp(A, perm, &sign); // LU decomp
        flag = gsl_linalg_LU_solve(A, perm , b, aux); // Solve system

        // Copy x_n , y_n into placeholders
        par.x[i] = gsl_vector_get(aux , gsl_permutation_get(perm , 0));
        par.y[i] = gsl_vector_get(aux , gsl_permutation_get(perm , 1));
    }

    // Deallocate system
    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(aux);
    gsl_permutation_free(perm);
}



void evaluateParametrization(Orbit orb , Param par , double smin , double smax) {
    int order = par.n;
    int nPoints = orb.n;
    double s;

    for (int i = 0 ; i < nPoints ; ++i) {
        s = smin + (smax - smin)*i/(nPoints-1.0);

        orb.x[i] = par.x[0]; orb.y[i] = par.y[0];
        for (int j = 1 ; j < order ; ++j) {
            orb.x[i] += par.x[j]*pow(s,j);
            orb.y[i] += par.y[j]*pow(s,j);
        }
    }
}