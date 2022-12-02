#include <R.h>
 #include <math.h>
 void g_Pex14_30xmdhwk ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[0]+p[1]*x[0+i**k] ;
y[1+i**l] = p[2]+p[1]*x[1+i**k] ;
y[2+i**l] = p[3]+p[1]*x[2+i**k] ;
y[3+i**l] = p[4]+p[1]*x[3+i**k] ;
y[4+i**l] = p[5]+p[1]*x[4+i**k] ;
y[5+i**l] = p[6]+p[1]*x[5+i**k] ;
y[6+i**l] = p[7]+p[1]*x[6+i**k] ;
y[7+i**l] = p[8]+p[1]*x[7+i**k] ;
y[8+i**l] = p[9]+p[1]*x[8+i**k] ; 
}
}