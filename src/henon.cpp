
#include "henon.h"

void f(Point* fp , Point p , double epsilon) {
    fp->x = p.x + epsilon*(p.y + epsilon*p.x*(1.0-p.x));
    fp->y = p.y + epsilon*p.x*(1.0-p.x);
}

void Df(gsl_matrix* M , Point p , double epsilon) {
    gsl_matrix_set(M , 0 , 0 , 1 + epsilon*epsilon*(1.0 - 2.0*p.x));
    gsl_matrix_set(M , 1 , 0 , epsilon);
    gsl_matrix_set(M , 0 , 1 , epsilon*(1.0 - 2.0*p.x));
    gsl_matrix_set(M , 1 , 1 , 1 );
}

void finv(Point* fp , Point p , double epsilon) {
    double aux = p.x - epsilon*p.y;
    fp->x = aux;
    fp->y = p.y - epsilon*aux*(1.0-aux);
}

void computeOrbitForward(Orbit orb , int n , Point p0 , double epsilon) {
    Point pk , pkp1;

    orb.n = n;
    orb.x[0] = p0.x; orb.y[0] = p0.y;
    pk.x = p0.x; pk.y = p0.y;

    for (int i = 1 ; i < n ; ++i){
        f(&pkp1 , pk , epsilon);
        orb.x[i] = pkp1.x; orb.y[i] = pkp1.y;
        pk.x = pkp1.x ; pk.y = pkp1.y;
    }
}

void computeOrbitBackward(Orbit orb , int n , Point p0 , double epsilon) {
    Point pk , pkp1;

    orb.n = n;
    orb.x[0] = p0.x; orb.y[0] = p0.y;
    pk.x = p0.x; pk.y = p0.y;

    for (int i = 1 ; i < n ; ++i){
        finv(&pkp1 , pk , epsilon);
        orb.x[i] = pkp1.x; orb.y[i] = pkp1.y;
        pk.x = pkp1.x ; pk.y = pkp1.y;
    }
}


