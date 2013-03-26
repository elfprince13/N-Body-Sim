#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef void (*TableauGen)(double *a, double *b, double *c);
typedef struct{ double *a; double *b; double *c; TableauGen init_tableau; } RKMethod;

void init_GLRK4(double *a, double *b, double *c){
	a[0][0] = 1/4.;
	a[0][1] = 1/4. - sqrt(3) / 6;
	a[1][0] = 1/4. + sqrt(3) / 6;
	a[1][1] = 1/4.;

	b[0] = 1/2.;
	b[1] = 1/2.;

	c[0] = a[0][0]+a[0][1];
	c[1] = a[1][0]+a[1][1];
}

double stepRungeKutta(double * yn,
