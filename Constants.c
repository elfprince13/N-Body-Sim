//
//  Constants.c
//  GravSim
//
//  Created by Thomas Dickerson on 4/15/13.
//
//

#include <stdio.h>
#include <stdlib.h>
#include "Constants.h"
#include <math.h>



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

double max(double f1, double f2){ return (f1 > f2)*f1 + (f1 <= f2)*f2;	}
double min(double f1, double f2){ return (f1 > f2)*f2 + (f1 <= f2)*f1;	}

void copyVector(double * data, Vector * v, int * index)
{	for(int i = 0; i < MAX_DIMS; (*index)++, i++){	data[(*index)] = getVField(v,i);	}	}
double smaller(double f){ return f * 1e-10; }
void copyVectorSmaller(double * data, Vector * v, int * index)
{	Vector sml = scale(smaller(1),*v);
	copyVector(data, &sml, index);	}

double ** makeDynamic2DArray(int r, int c){
	double ** parray = malloc(r*(sizeof(double *)+c*sizeof(double)));
	double *darray = (double*)(&(parray[r]));
	for(int i = 0; i < r; i++){
		parray[i] = &(darray[i*c]);
	}
	return parray;
}
