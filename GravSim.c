#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DIMS 2
#define MAX_DIMS 3


// Assume this is correct for now :)
#define EPSILON 2.22044604925031308e-16
#define sqEP 1.4901161193847656e-08
// in m^3 / (kg * s^2)
#define G 6.67384e-11

#define ALL_BODIES -1

#define X_START 0
#define P_START 3
#define MASS 6
#define RADIUS 7

typedef struct{double e1;double e2;double e3;} Vector;
typedef struct{Vector x; Vector p; double m; double r;} Body;
typedef double (*NHamiltonian)(Body ** bodies, int n, int which);
typedef struct{NHamiltonian hamiltonian; int n; Body ** bodies;} State;

double * getBFieldRef(Body * b,int field){
	double * ret;
	switch(field){
		case 0: ret = &(b->x.e1); break;
		case 1: ret = &(b->x.e2); break;
		case 2: ret = &(b->x.e3); break;
		case 3: ret = &(b->p.e1); break;
		case 4: ret = &(b->p.e2); break;
		case 5: ret = &(b->p.e3); break;
		case 6: ret = &(b->m); break;
		case 7: ret = &(b->r); break;
		default:
			printf("Error in getBFieldRef(): %d is not a valid field specifier\n",field);
			ret = NULL;
	}
	return ret;
}

double	getBField(Body * b,int field){	return *getBFieldRef(b, field);	}
void	setBField(Body * b,int field, double val){	*getBFieldRef(b, field) = val;	}
void	difBField(Body * b,int field, double val){	*getBFieldRef(b, field) += val;	}

double * getVFieldRef(Vector * v,int field){
	double * ret;
	switch(field){
		case 0: ret = &(v->e1); break;
		case 1: ret = &(v->e2); break;
		case 2: ret = &(v->e3); break;
		default:
			printf("Error in getVFieldRef(): %d is not a valid field specifier\n",field);
			ret = NULL;
	}
	return ret;
}

double	getVField(Vector * v,int field){	return *getVFieldRef(v, field);	}
void	setVField(Vector * v,int field, double val){	*getVFieldRef(v, field) = val;	}
void	difVField(Vector * v,int field, double val){	*getVFieldRef(v, field) += val;	}

Vector sum(const Vector v1, const Vector v2){	return (Vector){.e1=v1.e1+v2.e1,.e2=v1.e2+v2.e2,.e3=v1.e3+v2.e3};	}
Vector diff(const Vector v1, const Vector v2){	return (Vector){.e1=v1.e1-v2.e1,.e2=v1.e2-v2.e2,.e3=v1.e3-v2.e3};}
Vector scale(const double s, const Vector v1){	return (Vector){.e1=s*v1.e1,.e2=s*v1.e2,.e3=s*v1.e3};	}
double dot(const Vector v1, const Vector v2){	return (v1.e1*v2.e1) + (v1.e2*v2.e2) + (v1.e3*v2.e3);	}
double mag(const Vector v1){	return sqrt(v1.e1 * v1.e1 + v1.e2 * v1.e2 + v1.e3 * v1.e3);	}
Vector norm(const Vector v1){	return scale(1/mag(v1),v1);	}
const char * stringify(const Vector v1, char * buf){	sprintf(buf,"%G %G %G",v1.e1,v1.e2,v1.e3); return buf;	}

Vector gradient(State s,int n,int kind)
{
	Vector grad = (Vector){.e1=0,.e2=0,.e3=0};
	double H_l,H_h,h,e;
	Body * b = s.bodies[n];
	Body midpoint = *b;
	Body m_l, m_h;
	for(int i = kind; i < kind + DIMS; i++){
		h = fabs(sqEP * (sqEP + getBField(b,i)));		
		m_l = midpoint;
		m_h = midpoint;
		difBField(&m_l,i,-h);
		difBField(&m_h,i,h);
		s.bodies[n] = &m_l;
		
		H_l = s.hamiltonian(s.bodies,s.n,n);
		s.bodies[n] = &m_h;
		H_h = s.hamiltonian(s.bodies,s.n,n);
		
		setVField(&grad,i%MAX_DIMS,(H_h - H_l) / (2*h));
		s.bodies[n] = b;
	}
	return grad;
}

// The selector allows us to determine whether we want a single body or the whole system
// This should speed up differentiation by an order of N, since we may neglect terms
// which will not vary.
double newtonianGravitation(Body * bodies[], int n, int which){
	double H = 0;
	double H_i = 0;
	int istart = (which < 0) ? 0 : which;
	int iend = (which < 0) ? n : which + 1;
	int jend, i, j;
	for(i = istart; i < iend; i++){
		H_i = 0;
		jend = (which < 0) ? i : n;
		for(j  = !i ? 1 : 0; j < jend; j = (j+1 == i) ? j + 2 : j + 1){
			H_i += bodies[i]->m * bodies[j]->m / mag(diff(bodies[i]->x,bodies[j]->x));
		}
		H_i = dot(bodies[i]->p,bodies[i]->p) / (2 * bodies[i]->m) - G*(H_i);		
		H += H_i;
	}
	return H;
}


State buildSystem(double m[], double r[], Vector x[], Vector p[], int n){
	Body ** bodies = (Body**)calloc(sizeof(Body*),n);
	for(int i = 0; i < n; i++)	*(bodies[i] = (Body*)calloc(sizeof(Body),1)) = (Body){.x=x[i],.p=p[i],.m=m[i],.r=r[i]};
	return (State){.hamiltonian = &newtonianGravitation,.n = n,.bodies = bodies};
}

int teardownSystem(State s){
	for(int i = 0; i < s.n; i++) free(s.bodies[i]);
	free(s.bodies);
	s.bodies = NULL;
	return 0;
}

int main(void) {
	double mass[] = {1.9891e30,1.8986e27,3.3022e23};
	double radius[] = {696342e3,69911e3,2439.7e3};
	Vector position[] = {(Vector){.e1=0,.e2=0,.e3=0},(Vector){.e1=778547200e3,.e2=0,.e3=0},(Vector){.e1=57909100e3,.e2=0,.e3=0}};
	Vector momentum[] = {(Vector){.e1=0,.e2=0,.e3=0},(Vector){.e1=0,.e2=13.07e3 * mass[1],.e3=0},(Vector){.e1=0,.e2=47.87e3 * mass[2],.e3=0}};
	int body_count = 3;
	
	Vector xdot[body_count],pdot[body_count];
	State system = buildSystem(mass,radius,position,momentum,body_count);
	// Take some steps and see how it goes. 3.74336e8 seconds is Jupiter's orbit
	double dt = .1;
	double t = 0;
	double end_time = /*100000000;*/3.74336e8; //for jupiter or 7.6005e6 for mercury;
	double h_init = system.hamiltonian(system.bodies,system.n,ALL_BODIES);
	double h_drift, h_tot;
	int oline = 0;
	int operiod = (int)round(100000/dt);
	int i,j;
	char buf[6][64];
	for(t = 0; t < end_time; t += dt, oline = (oline + 1) % operiod){
		h_tot = system.hamiltonian(system.bodies,system.n,ALL_BODIES);
		h_drift = h_init - h_tot;
		if(!oline){
			printf("%G %G %G %s %s %s %s %s %s\n",t,h_tot,h_drift,
				   stringify(system.bodies[0]->x,buf[0]),stringify(system.bodies[0]->p,buf[1]),
				   stringify(system.bodies[1]->x,buf[2]),stringify(system.bodies[1]->p,buf[3]),
				   stringify(system.bodies[2]->x,buf[4]),stringify(system.bodies[2]->p,buf[5]));
		}
		for(i = 0; i < body_count; i++){
			xdot[i] = scale(dt,gradient(system,i,P_START));
			pdot[i] = scale(-dt,gradient(system,i,X_START));
		}
		for(i = 0; i < body_count; i++){
			for(j = 0; j < DIMS; j++){
				difBField(system.bodies[i],X_START+j,getVField(&xdot[i],j));
				difBField(system.bodies[i],P_START+j,getVField(&pdot[i],j));
			}
		}
	}
	return teardownSystem(system);
}
