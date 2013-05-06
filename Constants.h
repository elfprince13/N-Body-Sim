//
//  Constants.h
//  GravSim
//
//  Created by Thomas Dickerson on 4/15/13.
//
//

#ifndef GravSim_Constants_h
#define GravSim_Constants_h

#define DIMS 2
#define MAX_DIMS 3
#define TICK_MEMORY 1000/*000*/
#define TICKS_PER 100


// Assume this is correct for now :)
#define EPSILON 2.22044604925031308e-16
#define sqEP 1.4901161193847656e-08
// in m^3 / (kg * s^2)
#define G 6.67384e-11

#define ALL_BODIES -1

#define X_AND_P -1
#define X_START 0
#define P_START 3
#define MASS 6
#define RADIUS 7


#define LEFTH 0
#define RIGHTH 4
#define UPPERH 2
#define LOWERH 0
#define ALLQ 8

//--------------------------------------------------
// Typedefs
//--------------------------------------------------
typedef struct{double e1;double e2;double e3;} Vector;
typedef struct{Vector x; Vector p; double m; double r;} Body;

// Forward declaration
struct _state;
typedef Vector (*NGradient)(struct _state * s, int which, int kind);
typedef Body (*HamiltonsEquations)(struct _state * s, int which, Body * testpos, int kind);
typedef double (*NHamiltonian)(struct _state * s, int which);
typedef struct _state{NHamiltonian hamiltonian; NGradient analyticalGradient; int n; Body ** bodies;} State;

#define ZERO_VECTOR (Vector){.e1=0, .e2=0, .e3=0}



double * getBFieldRef(Body * b,int field);
double	getBField(Body * b,int field);
void	setBField(Body * b,int field, double val);
void	difBField(Body * b,int field, double val);
double * getVFieldRef(Vector * v,int field);

double	getVField(Vector * v,int field);
void	setVField(Vector * v,int field, double val);
void	difVField(Vector * v,int field, double val);

Vector sum(const Vector v1, const Vector v2);
Vector diff(const Vector v1, const Vector v2);
Vector scale(const double s, const Vector v1);
double dot(const Vector v1, const Vector v2);
double mag(const Vector v1);
Vector norm(const Vector v1);
const char * stringify(const Vector v1, char * buf);

double max(double f1, double f2);
double min(double f1, double f2);

void copyVector(double * data, Vector * v, int * index);
double smaller(double f);
void copyVectorSmaller(double * data, Vector * v, int * index);

double ** makeDynamic2DArray(int d1, int d2);

#endif
