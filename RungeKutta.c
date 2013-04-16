#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "RungeKutta.h"
#include "GravSim.h"

void free_GLRK4(void * tableau){
	NonSingularImplicitTableau * nsit = (NonSingularImplicitTableau *)tableau;
	if(nsit->a != NULL) free(nsit->a);
	if(nsit->b != NULL) free(nsit->b);
	if(nsit->c != NULL) free(nsit->c);
	if(nsit->d != NULL) free(nsit->d);
	free(nsit);
}

void free_StandardRK(void * tableau){
	NonSingularImplicitTableau * nsit = (NonSingularImplicitTableau *)tableau;
	if(nsit->a != NULL) free(nsit->a);
	if(nsit->b != NULL) free(nsit->b);
	if(nsit->c != NULL) free(nsit->c);
	free(nsit);
}
void free_PartitionedRK(void * tableau){
	NonSingularImplicitTableau * nsit = (NonSingularImplicitTableau *)tableau;
	if(nsit->a != NULL) free(nsit->a);
	if(nsit->b != NULL) free(nsit->b);
	if(nsit->c != NULL) free(nsit->c);
	free(nsit);
}

void init_GLRK4(void * tableau){
	double **a = ((NonSingularImplicitTableau *)tableau)->a;
	double *b = ((NonSingularImplicitTableau *)tableau)->b;
	double *c = ((NonSingularImplicitTableau *)tableau)->c;
	double *d = ((NonSingularImplicitTableau *)tableau)->d;
	
	a[0][0] = 1/4.;
	a[0][1] = 1/4. - sqrt(3) / 6;
	a[1][0] = 1/4. + sqrt(3) / 6;
	a[1][1] = 1/4.;

	b[0] = 1/2.;
	b[1] = 1/2.;

	c[0] = a[0][0]+a[0][1];
	c[1] = a[1][0]+a[1][1];
	
	double inv[2][2] = {{0,0},{0,0}};
	double det = (a[0][0]*a[1][1])-(a[0][1]*a[1][0]);
	if( fabs(det) < 2 * EPSILON){
		printf("Error, singluar A matrix in init_GLRK4\n");
	} else{
		inv[0][0]=a[1][1]/det;
		inv[1][1]=a[0][0]/det;
		inv[0][1]=-a[0][1]/det;
		inv[1][0]=-a[1][0]/det;
		d[0]=(b[0]*inv[0][0]+b[1]*inv[1][0]);
		d[1]=(b[0]*inv[0][1]+b[1]*inv[0][1]);
	}
}

void init_ForwardEuler(void * tableau){
	double **a = ((StandardTableau *)tableau)->a;
	double *b = ((StandardTableau *)tableau)->b;
	double *c = ((StandardTableau *)tableau)->c;
	a[0][0] = 0;
	
	b[0] = 1;
	
	c[0] = a[0][0];
	
}

void init_RK4(void * tableau){
	double **a = ((StandardTableau *)tableau)->a;
	double *b = ((StandardTableau *)tableau)->b;
	double *c = ((StandardTableau *)tableau)->c;

	a[0][0] = 0;	a[0][1] = 0;	a[0][2] = 0;	a[0][3] = 0;
	a[1][0] = 1/2.;	a[1][1] = 0;	a[1][2] = 0;	a[1][3] = 0;
	a[2][0] = 0;	a[2][1] = 1/2.;	a[2][2] = 0;	a[2][3] = 0;
	a[3][0] = 0;	a[3][1] = 0;	a[3][2] = 1;	a[3][3] = 0;
	
	b[0] = 1/6.;
	b[1] = 1/3.;
	b[2] = 1/3.;
	b[3] = 1/6.;
	
	c[0] = a[0][0]+a[0][1]+a[0][2]+a[0][3];
	c[1] = a[1][0]+a[1][1]+a[1][2]+a[1][3];
	c[2] = a[2][0]+a[2][1]+a[2][2]+a[2][3];
	c[3] = a[3][0]+a[3][1]+a[3][2]+a[3][3];
}

void init_PRK6(void * tableau){
	double **a = ((PartitionedTableau *)tableau)->a;
	double *b = ((PartitionedTableau *)tableau)->b;
	double *c = ((PartitionedTableau *)tableau)->c;
	double **A = ((PartitionedTableau *)tableau)->A;
	double *B = ((PartitionedTableau *)tableau)->B;
	double *C = ((PartitionedTableau *)tableau)->C;
	
	double d_[3];
	double c_[3];
	
	// +real solution to 12*z^4 - 24*z^2 + 16*z - 3 = 0
	/*d[0] = (1/2.) * sqrt(8/3.+1/(3.*pow((2-sqrt(3)),(1/3.)))+ 1/3. * pow((2-sqrt(3)),(1/3.)) +8/sqrt(3 * (4-1/pow((2-sqrt(3)),(1/3.)) -pow((2-sqrt(3)),(1/3.)))))	-1/(2. * sqrt(3/(4-1/pow((2-sqrt(3)),(1/3.))-pow((2-sqrt(3)),(1/3.))))); */
	c_[2] = d_[0] = 0.9196615230173998570508976381533827895633;
	c_[1] = d_[1] = 1.5050911349195523119693694258934572;
	c_[0] = d_[2] = 1 - d_[0] - d_[1];
	
	b[0] = c_[0]/2;
	b[1] = c_[1]/2;
	b[2] = c_[2]/2;
	b[3] = c_[2]/2;
	b[4] = c_[1]/2;
	b[5] = c_[0]/2;
	
	B[0] = d_[0]/2;
	B[1] = d_[1]/2;
	B[2] = d_[2]/2;
	B[3] = d_[2]/2;
	B[4] = d_[1]/2;
	B[5] = d_[0]/2;
	
	for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
			if( i < 3){
				a[i][j] = (j >= i) ? b[i] : 0;
				A[i][j] = (j > i) ? B[i] : 0;
			} else{
				a[i][j] = (j > i) ? b[i] : 0;
				A[i][j] = (j >= i) ? B[i] : 0;
			}
			c[j] += a[i][j];
			C[j] += A[i][j];
		}
		
	}
}

Body scaleDeltaBody(double s,Body * b){
	Body ret = *b;
	ret.x = scale(s,ret.x);
	ret.p = scale(s,ret.p);
	return ret;
}

Body integrate_RK(State * s, int which, ODE_RHS f, RKMethod * rkm){
	int stages = rkm->stages;
	Body k[stages],passthrough,*ibp,ib,tmp_delta;
	ibp = s->bodies[which];
	ib = *ibp;
	Body ret = (Body){.x = ZERO_VECTOR, .p = ZERO_VECTOR, .m = ib.m, .r = ib.r};
	
	StandardTableau * tableau = (StandardTableau *)(rkm->tableau);
	
	double **a = (tableau)->a;
	double *b = (tableau)->b;
	double *c = (tableau)->c;
	
	// Autonomous, so we can neglect c!
	for(int i = 0; i < stages; i++){
		passthrough = ib;
		for(int j = 0; j < stages; j++){
			//printf("%f ",a[i][j]);
			tmp_delta = scaleDeltaBody(dt*a[i][j],ibp);
			passthrough.x = sum(passthrough.x,tmp_delta.x);
			passthrough.p = sum(passthrough.p,tmp_delta.p);
		}
		//printf("\n");
		k[i] = f(s, which, &passthrough);//&(scaleDeltaBody(&k[i],b[i]));
		passthrough = scaleDeltaBody(dt*b[i],&(k[i]));
		ret.x = sum(ret.x,passthrough.x);
		ret.p = sum(ret.p,passthrough.p);
	}
	//printf("\n");
	return ret;
		
}

Body integrate_PRK(State * s, int which, ODE_RHS f, RKMethod * rkm){
	int stages = rkm->stages;
	Body k[stages],klc[stages],kuc[stages],passthrough,passthroughlc,passthroughuc,*ibp,ib;
	Vector tmp_deltalc,tmp_deltauc;
	ibp = s->bodies[which];
	ib = *ibp;
	Body ret = (Body){.x = ZERO_VECTOR, .p = ZERO_VECTOR, .m = ib.m, .r = ib.r};
	
	PartitionedTableau * tableau = (PartitionedTableau *)(rkm->tableau);
	
	double **a = (tableau)->a;
	double *b = (tableau)->b;
	double *c = (tableau)->c;
	double **A = (tableau)->A;
	double *B = (tableau)->B;
	double *C = (tableau)->C;
	
	// Autonomous, so we can neglect c!
	for(int i = 0; i < stages; i++){
		passthrough = ib;
		for(int j = 0; j < stages; j++){
			tmp_deltauc = scale(dt*A[i][j],ib.x);//scaleDeltaBody(dt*a[i][j],ibp);
			tmp_deltalc = scale(dt*a[i][j],ib.p);
			passthrough.x = sum(passthrough.x,tmp_deltauc);
			passthrough.p = sum(passthrough.p,tmp_deltalc);
		}
		k[i] = f(s, which, &passthrough);//&(scaleDeltaBody(&k[i],b[i]));
		passthrough.x = scale(dt*B[i],(k[i].x));
		passthrough.p = scale(dt*b[i],(k[i].p));
		ret.x = sum(ret.x,passthrough.x);
		ret.p = sum(ret.p,passthrough.p);
	}
	return ret;
	
}
