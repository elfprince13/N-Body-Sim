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
	StandardTableau * st = (StandardTableau *)tableau;
	if(st->a != NULL) free(st->a);
	if(st->b != NULL) free(st->b);
	if(st->c != NULL) free(st->c);
	free(st);
}
void free_PartitionedRK(void * tableau){
	PartitionedTableau * pt = (PartitionedTableau *)tableau;
	if(pt->a != NULL) free(pt->a);
	if(pt->b != NULL) free(pt->b);
	/*
	if(pt->c != NULL) free(pt->c);
	if(pt->A != NULL) free(pt->A);
	if(pt->B != NULL) free(pt->B);
	if(pt->C != NULL) free(pt->C);
	 //*/
	free(pt);
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
	/*
	double **a = ((PartitionedTableau *)tableau)->a;
	double *b = ((PartitionedTableau *)tableau)->b;
	double *c = ((PartitionedTableau *)tableau)->c;
	double **A = ((PartitionedTableau *)tableau)->A;
	double *B = ((PartitionedTableau *)tableau)->B;
	double *C = ((PartitionedTableau *)tableau)->C;
	
	double d_[3];
	double c_[3];
	 //*/
	// +real solution to 12*z^4 - 24*z^2 + 16*z - 3 = 0
	 /*d[0] = (1/2.) * sqrt(8/3.+1/(3.*pow((2-sqrt(3)),(1/3.)))+ 1/3. * pow((2-sqrt(3)),(1/3.)) +8/sqrt(3 * (4-1/pow((2-sqrt(3)),(1/3.)) -pow((2-sqrt(3)),(1/3.)))))	-1/(2. * sqrt(3/(4-1/pow((2-sqrt(3)),(1/3.))-pow((2-sqrt(3)),(1/3.))))); */
	/*
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
				a[i][j] = (i >= j) ? b[j] : 0;
				A[i][j] = (i > j) ? B[j] : 0;
			} else{
				a[i][j] = (i > j) ? b[j] : 0;
				A[i][j] = (i >= j) ? B[j] : 0;
			}
			c[j] += a[i][j];
			C[j] += A[i][j];
			printf("%.3f (%.3f)\t",a[i][j],A[i][j]);
		}
		printf("\n");
		
		
	}
	//*/
	double *a = ((PartitionedTableau *)tableau)->a;
	double *b = ((PartitionedTableau *)tableau)->b;
	
	/*
	a[0] = 0.0502627644003922;
	a[1] = 0.413514300428344;
	a[2] = 0.0450798897943977;
	a[3] = -0.188054853819569;
	a[4] = 0.541960678450780;
	a[5] = 1 - 2*(a[0]+a[1]+a[2]+a[3]+a[4]);

	
	b[0] = 0.148816447901042;
	b[1] = -0.132385865767784;
	b[2] = 0.067307604692185;
	b[3] = 0.432666402578175;
	b[4] = 0.5 - (b[0]+b[1]+b[2]+b[3]);
	b[5] = 0;
	printf("%f %f\n",a[0]+a[1]+a[2]+a[3]+a[4]+a[5],b[0]+b[1]+b[2]+b[3]+b[4]+b[5]);
	//*/
	
	// Ruth 1990
	a[0] = 0.339839625839110000;
	a[1] =-0.088601336903027329;
	a[2] = 0.5858564768259621188;
	a[3] =-0.603039356536491888;
	a[4] = 0.3235807965546976394;
	a[5] = 0.4423637942197494587;
	
	b[0] = 0.1193900292875672758;
	b[1] = 0.6989273703824752308;
	b[2] =-0.1713123582716007754;
	b[3] = 0.4012695022513534480;
	b[4] = 0.0107050818482359840;
	b[6] =-0.0589796254980311632;
	//printf("hi\n");
	//double as = a[0]+a[1]+a[2]+a[3]+a[4]+a[5];
	//double bs = b[0]+b[1]+b[2]+b[3]+b[4]+b[5];
	//printf("%f %f\n",0.f,0.f);
	
	
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
		//printf("k_%d: ",i);
		for(int j = 0; j < i; j++){
			//printf("%f ",a[i][j]);
			tmp_delta = scaleDeltaBody(dt*a[i][j],&(k[j])); //ibp);
			passthrough.x = sum(passthrough.x,tmp_delta.x);
			passthrough.p = sum(passthrough.p,tmp_delta.p);
		}
		//printf("\n");
		k[i] = f(s, which, &passthrough,X_AND_P);//&(scaleDeltaBody(&k[i],b[i]));
		passthrough = scaleDeltaBody(dt*b[i],&(k[i]));
		ret.x = sum(ret.x,passthrough.x);
		ret.p = sum(ret.p,passthrough.p);
	}
	//printf("\n");
	return ret;
		
}

/*
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
		//passthroughuc = ib;
		//passthroughlc = ib;
		passthrough = ib;
		for(int j = 0; j <= i; j++){
			//tmp_deltauc = scale(dt*A[i][j],(kuc[j].x));//ib.x);//scaleDeltaBody(dt*a[i][j],ibp);
			//tmp_deltalc = scale(dt*a[i][j],(klc[j].p));//ib.p);
			tmp_deltauc = scale(dt*A[i][j],(k[j].x));//ib.x);//scaleDeltaBody(dt*a[i][j],ibp);
			tmp_deltalc = scale(dt*a[i][j],(k[j].p));//ib.p);
			//passthroughuc.x = sum(passthroughuc.x,tmp_deltauc);
			//passthroughlc.p = sum(passthroughlc.p,tmp_deltalc);
			passthrough.x = sum(passthrough.x,tmp_deltauc);
			passthrough.p = sum(passthrough.p,tmp_deltalc);
		}
		//kuc[i] = f(s, which, &passthroughuc);//&(scaleDeltaBody(&k[i],b[i]));
		//klc[i] = f(s, which, &passthroughlc);
		k[i] = f(s, which, &passthrough,X_AND_P);
		//passthrough.x = scale(dt*B[i],(kuc[i].x));
		//passthrough.p = scale(dt*b[i],(klc[i].p));
		passthrough.x = scale(dt*B[i],(k[i].x));
		passthrough.p = scale(dt*b[i],(k[i].p));
		ret.x = sum(ret.x,passthrough.x);
		ret.p = sum(ret.p,passthrough.p);
	}
	return ret;
	
}
 //*/


Body integrate_PRK(State * s, int which, ODE_RHS f, RKMethod * rkm){
	int stages = rkm->stages;
	Body k[stages+2],ib,*ibp,tmp;
	ibp = s->bodies[which];
	ib = *ibp;
	Body ret = (Body){.x = ZERO_VECTOR, .p = ZERO_VECTOR, .m = ib.m, .r = ib.r};
	
	PartitionedTableau * tableau = (PartitionedTableau *)(rkm->tableau);
	
	double *a = (tableau)->a;
	double *b = (tableau)->b;
	double as = a[0]+a[1]+a[2]+a[3]+a[4]+a[5];
	double bs = b[0]+b[1]+b[2]+b[3]+b[4]+b[5];
	
	k[0] = ib;
	k[1] = ib;
	for(int i = 1; i <= stages; i++){
		// This should be safe, as f(,,,X_START) depends only on p[i-1]
		// and f(,,,P_START) depends only on x[i]
		k[i+1] = k[i];
		
		tmp = f(s,which,&(k[i]), P_START);
		tmp = scaleDeltaBody(dt*a[i-1], &tmp);
		k[i].p = sum(k[i-1].p,tmp.p);
		
		tmp = f(s,which,&(k[i]), X_START);
		tmp = scaleDeltaBody(dt*b[i-1], &tmp);
		k[i+1].x = sum(k[i].x,tmp.x);
	}
	//return k[stages];
	ret.p = diff(k[stages].p,ib.p);
	ret.x = diff(k[stages+1].x,ib.x);
	return ret;
	
}
