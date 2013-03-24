#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DIMS 2
#define MAX_DIMS 3


// Assume this is correct for now :)
#define EPSILON 2.22044604925031308e-16
#define sqEP 1.4901161193847656e-08

typedef struct{double e1;double e2;double e3;} Vector;

double * getVFieldRef(Vector * v,int field){
	double * ret;
	switch(field){
		case 0:
			ret = &(v->e1); break;
		case 1:
			ret = &(v->e2); break;
		case 2:
			ret = &(v->e3); break;
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


typedef double (*Surface)(Vector pos);

double vsin(Vector pos)
{
	return sin(mag(pos)); 
}

double invsq(Vector pos)
{
	return 1/(mag(pos)*mag(pos));
}

Vector ngrad(Surface s, Vector x)
{
	Vector grad = (Vector){.e1=0,.e2=0,.e3=0};
	int i;
	double S_l,S_h,h,e;
	Vector x_l, x_h;
	char buf[64];
	for(i = 0; i < DIMS; i++){
		h = fabs(sqEP * (sqEP + getVField(&x,i)));		
		x_l = x;
		x_h = x;
		difVField(&x_l,i,-h);
		difVField(&x_h,i,h);
		
		S_l = s(x_l);
		S_h = s(x_h);
		
		setVField(&grad,i,(S_h - S_l) / (2*h));
	}
	return grad;
}

int main(void) {
	int surfcount = 3;
	Surface surfs[] = {&vsin,&mag,&invsq};
	double i,j;
	int k;
	char buf[64];
	for(i=-10;i<10.5;i+=.5) {
		for(j=-10;j<10.5;j+=.5) {
			Vector p = (Vector){.e1=i,.e2=j,.e3=0};
			printf("%s", stringify(p,buf));
			for(k=0;k<surfcount; k++){
				Vector g = ngrad(surfs[k],p);
				printf(" %s", stringify(scale(0.5,g),buf));
			}
			printf("\n");
		}
	}
	return 0;
}