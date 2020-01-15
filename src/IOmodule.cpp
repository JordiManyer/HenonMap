

#include "henon.h"

// Print in terminal
void printOrbit(Orbit orb){
    for (int i=0;i<orb.n;++i){
        printf("%g %g \n",orb.x[i],orb.y[i]);
    }
}

void printPoint(Point p){
    printf("%g %g \n",p.x,p.y);
}

void printArray(int n , double* v){
    for (int i = 0 ; i < n ; ++i) printf("%g, " , v[i]);
    printf("\n");
}



// Print to file
void tofileOrbit(FILE * fptr, Orbit orb){
    fprintf(fptr , "%d\n", orb.n);
    for (int i=0;i<orb.n;++i){
        fprintf(fptr , "%12.6g %12.6g \n",orb.x[i],orb.y[i]);
    }
}