#include <iostream>
#include "src/henon.h"

int main() {
    int dim; double eps;
    SStab eq; Point p;
    gsl_complex lambda;
    FILE* fptr;

    // Parameters
    eps = 0.5;
    p.x = 0.0; p.y = 0.0;
    printf("Parameters : \n");
    printf("epsilon = %g , p0 = (%g , %g) \n\n", eps , p.x , p.y);

    // Allocate space for DF(0,0), eigenvalues and eigenvectors
    dim = 2;
    eq.n = dim; eq.M = gsl_matrix_calloc(dim,dim);
    eq.evec = gsl_matrix_complex_calloc(dim,dim); eq.eval = gsl_vector_complex_calloc(dim);
    int aux [dim]; eq.stable = aux;

    // Calculate DF(0,0) and it's evals, evects
    Df(eq.M , p , eps);
    getStability(eq);
    for (int i = 0 ; i < eq.n ; ++i){
        lambda = gsl_vector_complex_get(eq.eval , i);
        printf("Eigenvalue number %d : %g + i%g \n" , i , GSL_REAL(lambda) , GSL_IMAG(lambda));
    }

    // Analytical eigenvalues
    double ev1, ev2;
    ev1 = (2.0 + eps*eps + eps*sqrt(eps*eps+4))/2.0;
    ev2 = (2.0 + eps*eps - eps*sqrt(eps*eps+4))/2.0;
    printf("Analytical eigenvalues: %g %g \n" , ev1 , ev2);

    // Test eigenvectors
    gsl_matrix_complex* test;
    test = gsl_matrix_complex_calloc(eq.n , eq.n);
    for (int i = 0; i < eq.n; ++i) {
        for (int j = 0; j < eq.n; ++j) {
            for (int k = 0; k<eq.n; ++k){
                lambda = gsl_complex_rect(gsl_matrix_get(eq.M,i,k),0);
                lambda = gsl_complex_mul(lambda , gsl_matrix_complex_get(eq.evec,k,j));
                lambda = gsl_complex_add(lambda , gsl_matrix_complex_get(test,i,j));
                gsl_matrix_complex_set(test,i,j,lambda);
            }
        }
    }
    double error = 0.0;
    for (int i = 0; i < eq.n; ++i) {
        for (int j = 0; j < eq.n; ++j) {
            lambda = gsl_matrix_complex_get(test,i,j);
            lambda = gsl_complex_div(lambda , gsl_vector_complex_get(eq.eval,j));
            gsl_matrix_complex_set(test,i,j,lambda);
            error += gsl_complex_abs(gsl_complex_sub(lambda,gsl_matrix_complex_get(eq.evec,i,j)));
        }
    }
    printf("Are eigenvectors correct? %d , Error = %g \n" , gsl_matrix_complex_equal(test , eq.evec) , error);


    // Get parametrizations
    Param parS , parU;
    int order = 200;
    parS.n = order; parU.n = order;
    double xs[order], ys[order]; parS.x = xs; parS.y = ys;
    double xu[order], yu[order]; parU.x = xu; parU.y = yu;

    getParametrization(parS , parU , dim , p , eq , eps);

    printf("\n Stable parametrization: \n ");
    for (int i = 0; i < order; ++i) printf("%g " , parS.x[i]);
    printf("\n");
    for (int i = 0; i < order; ++i) printf("%g " , parS.y[i]);
    printf("\n Unstable parametrization: \n ");
    for (int i = 0; i < order; ++i) printf("%g " , parU.x[i]);
    printf("\n");
    for (int i = 0; i < order; ++i) printf("%g " , parU.y[i]);
    printf("\n");


    // Stable Manifold
    int nPoints = 50; double smin = 1.e-12; double smax = 1.0; // 1.e-7;
    Orbit orbS; orbS.n = nPoints;
    double xos[nPoints] , yos[nPoints]; orbS.x = xos; orbS.y = yos;
    evaluateParametrization(orbS , parS , smin , smax);

    fptr = fopen("Stable.out","w");
    tofileOrbit(fptr , orbS);
    fclose(fptr);

    fptr = fopen("Stable_orbits.out","w");
    fprintf(fptr , "%d \n" , nPoints);
    Point paux; Orbit orbaux; orbaux.n = nPoints;
    double xaux[nPoints] , yaux[nPoints]; orbaux.x = xaux; orbaux.y = yaux;
    for (int i = 0 ; i < nPoints ; ++i) {
        paux.x = orbS.x[i]; paux.y = orbS.y[i];
        computeOrbitBackward(orbaux , nPoints , paux , eps);
        tofileOrbit(fptr , orbaux);
    }
    fclose(fptr);


    // Unstable Manifold
    Orbit orbU; orbU.n = nPoints;
    double xou[nPoints] , you[nPoints]; orbU.x = xou; orbU.y = you;
    evaluateParametrization(orbU , parU , smin , smax);

    fptr = fopen("Unstable.out","w");
    tofileOrbit(fptr , orbU);
    fclose(fptr);

    fptr = fopen("Unstable_orbits.out","w");
    fprintf(fptr , "%d \n" , nPoints);
    orbaux.n = nPoints;
    for (int i = 0 ; i < nPoints ; ++i) {
        paux.x = orbU.x[i]; paux.y = orbU.y[i];
        computeOrbitForward(orbaux , nPoints , paux , eps);
        tofileOrbit(fptr , orbaux);
    }
    fclose(fptr);
}
