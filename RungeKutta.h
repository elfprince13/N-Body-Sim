//
//  RungeKutta.h
//  GravSim
//
//  Created by Thomas Dickerson on 4/15/13.
//
//

#ifndef GravSim_RungeKutta_h
#define GravSim_RungeKutta_h
#include <math.h>
#include "Constants.h"

#define ODE_RHS HamiltonsEquations
typedef void (*TableauCleaner)(void*);
typedef struct{ double **a; double *b; double *c; double *d; TableauCleaner tfree;} NonSingularImplicitTableau;
typedef struct{ double /***/*a; double *b; /*double *c; double **A; double *B; double *C;*/ TableauCleaner tfree;} PartitionedTableau;
typedef struct{ double **a; double *b; double *c; TableauCleaner tfree;} StandardTableau;

struct _rkm;

typedef void (*TableauGen)(void * tableau);
typedef Body (*Integrator)(State * s, int which, ODE_RHS f, struct _rkm * rkm);
typedef struct _rkm{ void * tableau; TableauGen init_tableau; Integrator integrate; int stages;} RKMethod;


void init_ForwardEuler(void *);
void init_GLRK4(void *);
void init_RK4(void *);
void init_PRK6(void *);
Body integrate_RK(State * s, int which, ODE_RHS f, RKMethod * rkm);
Body integrate_PRK(State * s, int which, ODE_RHS f, RKMethod * rkm);

void free_GLRK4(void * tableau);
void free_StandardRK(void * tableau);
void free_PartitionedRK(void * tableau);

#endif
