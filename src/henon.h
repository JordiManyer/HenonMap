
#ifndef PROJECTQQMDS_HENONMAP_HENON_H
#define PROJECTQQMDS_HENONMAP_HENON_H



#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


/////////////////////////////////////////
/*             STRUCTURES              */
/////////////////////////////////////////

struct parametrization {
    int n;
    double* x;
    double* y;
};
typedef struct parametrization Param;

struct point {
    double x;
    double y;
};
typedef struct point Point;

struct orbit {
    int n;
    double * x;
    double * y;
};
typedef struct orbit Orbit;

struct stabilityStruct{ // Stability features
    int n; // dimension of the matrix.
    gsl_matrix * M; // Matrix
    gsl_vector_complex * eval; // Vector of eigenvalues
    gsl_matrix_complex * evec; // Matrix which columns are eigenvectors
    int* stable; // Classification vector: 1->stable eigenvalue, 2->Unstable eigenvalue, 3->Pure complex eigenvalue
};
typedef struct stabilityStruct SStab;



/////////////////////////////////////////
/*      FUNCTIONS AND SUBROUTINES      */
/////////////////////////////////////////

// henon.cpp -> Definition of the map
void f(Point* fp , Point p , double epsilon);
void finv(Point* fp , Point p , double epsilon);
void Df(gsl_matrix* M , Point p , double epsilon);
void computeOrbitForward(Orbit orb , int n , Point p0 , double epsilon);
void computeOrbitBackward(Orbit orb , int n , Point p0 , double epsilon);

// stability.cpp -> Stability of eq points
void getStability(SStab M);

// paramMethod.cpp
void getParametrization(Param parStable , Param parUnst , int dim , Point x0 , SStab M , double eps);
void paramMethod1D(Param par , int dim , Point p0 , double eval , gsl_vector* evec , double eps);
void paramMethod1D_finv(Param par , int dim , Point p0 , double eval , gsl_vector* evec , double eps);
void evaluateParametrization(Orbit orb , Param par , double smin , double smax);

// IOmodule.cpp
void tofileOrbit(FILE * fptr, Orbit orb);
void printOrbit(Orbit orb);
void printPoint(Point p);
void printArray(int n , double* v);



#endif //PROJECTQQMDS_HENONMAP_HENON_H
