#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "zpr.h"

#define DIMS 2
#define MAX_DIMS 3
#define TICK_MEMORY 1000000
#define TICKS_PER 1000


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


static double _matrix[16];
GLdouble * getMatrix()
{
	glGetDoublev(GL_MODELVIEW_MATRIX,_matrix);
	return _matrix;
}

double determinant(const GLdouble *m)
{
/* NB. OpenGL Matrices are COLUMN major. */
#define MAT(m,r,c) (m)[(c)*4+(r)]

/* Here's some shorthand converting standard (row,column) to index. */
#define m11 MAT(m,0,0)
#define m12 MAT(m,0,1)
#define m13 MAT(m,0,2)
#define m14 MAT(m,0,3)
#define m21 MAT(m,1,0)
#define m22 MAT(m,1,1)
#define m23 MAT(m,1,2)
#define m24 MAT(m,1,3)
#define m31 MAT(m,2,0)
#define m32 MAT(m,2,1)
#define m33 MAT(m,2,2)
#define m34 MAT(m,2,3)
#define m41 MAT(m,3,0)
#define m42 MAT(m,3,1)
#define m43 MAT(m,3,2)
#define m44 MAT(m,3,3)

   GLdouble d12, d13, d23, d24, d34, d41;
   GLdouble tmp[16]; /* Allow out == in. */

   /* Inverse = adjoint / det. (See linear algebra texts.)*/

   /* pre-compute 2x2 dets for last two rows when computing */
   /* cofactors of first two rows. */
   d12 = (m31*m42-m41*m32);
   d13 = (m31*m43-m41*m33);
   d23 = (m32*m43-m42*m33);
   d24 = (m32*m44-m42*m34);
   d34 = (m33*m44-m43*m34);
   d41 = (m34*m41-m44*m31);

   tmp[0] =  (m22 * d34 - m23 * d24 + m24 * d23);
   tmp[1] = -(m21 * d34 + m23 * d41 + m24 * d13);
   tmp[2] =  (m21 * d24 + m22 * d41 + m24 * d12);
   tmp[3] = -(m21 * d23 - m22 * d13 + m23 * d12);

   /* Compute determinant as early as possible using these cofactors. */
   return m11 * tmp[0] + m12 * tmp[1] + m13 * tmp[2] + m14 * tmp[3];
}

//--------------------------------------------------
// Typedefs
//--------------------------------------------------
typedef struct{double e1;double e2;double e3;} Vector;
typedef struct{Vector x; Vector p; double m; double r;} Body;
typedef double (*NHamiltonian)(Body ** bodies, int n, int which);
typedef struct{NHamiltonian hamiltonian; int n; Body ** bodies;} State;

//--------------------------------------------------
// System parameters
//--------------------------------------------------
static State gravsys;
static int body_count = 3;

static double dt = 10;//.1;
static double t = 0;
static double end_time = /*100000000;*/ 3.74336e8; //for jupiter or 7.6005e6 for mercury;

static int operiod;
static int  oline;
static int ooperiod;
static int ooline;

static double h_init, h_drift;

static double * positions;
static double * sizes;
static int * bodyIndices;



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

double sgnLog(double f){ return copysign(max(log10(fabs(f)+EPSILON),0),f); }
void copyVectorLog(double *data, Vector * v, int *index){
	int i;
	for(i = 0; i < MAX_DIMS; (*index)++, i++){
		data[(*index)] = sgnLog(getVField(v,i));
	}
	//printf("%G %G -> %G %G @ %d\n",getVField(v,i-2),getVField(v,i-1),data[(*index)-2],data[(*index)-1],(*index)-2);
}

double sgnSqrt(double f){	return copysign(sqrt(fabs(f)),f);	}
void copyVectorSqrt(double * data, Vector * v, int * index){
	int i;
	for(i = 0; i < MAX_DIMS; (*index)++, i++){
		data[(*index)] = sgnSqrt(getVField(v,i));
	}
	//printf("%G %G -> %G %G @ %d\n",getVField(v,i-2),getVField(v,i-1),data[(*index)-2],data[(*index)-1],(*index)-2);
}

double smaller(double f){ return f * 1e-10; }
void copyVectorSmaller(double * data, Vector * v, int * index){
	int i;
	for(i = 0; i < MAX_DIMS; (*index)++, i++){
		data[(*index)] = smaller(getVField(v,i));
	}
	//printf("%G %G %G -> %G %G %G @ %d\n",getVField(v,i-3),getVField(v,i-2),getVField(v,i-1),data[(*index)-3],data[(*index)-2],data[(*index)-1],(*index)-2);
}


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

void stepSys(int n){
	double h_tot;
	Vector xdot[body_count],pdot[body_count];
	char buf[64][64];
	for(;t < end_time && n > 0; t += dt, oline = (oline + 1) % operiod, ooline = (ooline + !oline) % ooperiod, n--){
		h_tot = gravsys.hamiltonian(gravsys.bodies,gravsys.n,ALL_BODIES);
		h_drift = h_init - h_tot;
		if(!oline){
			/*printf("%G %G %G %s %s %s %s %s %s\n",t,h_tot,h_drift,
						stringify(gravsys.bodies[0]->x,buf[0]),stringify(gravsys.bodies[0]->p,buf[1]),
						stringify(gravsys.bodies[1]->x,buf[2]),stringify(gravsys.bodies[1]->p,buf[3]),
						stringify(gravsys.bodies[2]->x,buf[4]),stringify(gravsys.bodies[2]->p,buf[5]));*/
			
			if(!ooline){
				for(int i = 0; i < body_count; i++){
					bodyIndices[2*i] = bodyIndices[2*i+1];
				}	
			}
			
			for(int i = 0; i < body_count; i++){
				copyVectorSmaller(positions,&(gravsys.bodies[i]->x), &(bodyIndices[2*i]));
			}
		}
		for(int i = 0; i < body_count; i++){
			xdot[i] = scale(dt,gradient(gravsys,i,P_START));
			pdot[i] = scale(-dt,gradient(gravsys,i,X_START));
		}
		for(int i = 0; i < body_count; i++){
			for(int j = 0; j < DIMS; j++){
				difBField(gravsys.bodies[i],X_START+j,getVField(&xdot[i],j));
				difBField(gravsys.bodies[i],P_START+j,getVField(&pdot[i],j));
			}
		}
	}
}

#define GLERROR                                                    \
    {                                                              \
        GLenum code = glGetError();                                \
        while (code!=GL_NO_ERROR)                                  \
        {                                                          \
            printf("%s\n",(char *) gluErrorString(code));          \
                code = glGetError();                               \
        }                                                          \
    }

void drawSys()
{
	static float colors[3] = {0.9,0.3,0.6};
	double scale = pow(determinant(getMatrix()),1/3.);
	int toDraw = ((t / dt) / operiod > ooperiod) ? ooperiod : ooline;
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_CULL_FACE);
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(MAX_DIMS,GL_DOUBLE,0,positions);
	
	for(int i = 0; i < body_count; i++)
	{
		glPushName(i);
		glColor3f(colors[i%3],colors[(i+1)%3],colors[(i+2)%3]);
		glPointSize(sizes[i]*scale);
		glDrawArrays(GL_POINTS,i*TICK_MEMORY,toDraw);
		glPopName();
	}
	
	glDisableClientState(GL_VERTEX_ARRAY);

	
	glEnable(GL_CULL_FACE);
	glDisable(GL_BLEND);
	glDisable(GL_POINT_SMOOTH);
}

void drawTime()
{
	static char buf[128];
	sprintf(buf,"Time: %G / %G by %G",t,end_time,dt);
	glRasterPos2i(-100,100);
	glColor3f(1.0f, 1.0f, 0.3f);
	for(char *b=buf;*b;b++) glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *b);
}

/* Callback functions for drawing */

void tdisplay(int done)
{
   GLERROR;

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// printf("Stepping time from %f\n",t);
	stepSys(TICKS_PER);
	//printf("Drawing system\n");
	drawSys();
	drawTime();
	//printf("next frame\n\n");
   glutSwapBuffers();

   GLERROR;
	if(!done) glutTimerFunc(17,tdisplay,t >= end_time);
}
void display(){ tdisplay(1); /* Don't schedule extra timers */ }

/* Callback function for pick-event handling from ZPR */

void pick(GLint name)
{
   printf("Pick: %d\n",name);
   fflush(stdout);
}

/* Entry point */

int main(int argc, char *argv[])
{
	

	
    /* Initialise GLUT and create a window */

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800,600);
    glutCreateWindow("GravSim");

    /* Configure GLUT callback functions */

    glutDisplayFunc(display);
	glutTimerFunc(17,tdisplay,t >= end_time);

    glScalef(.008,.008,.008);

    /* Configure ZPR module */

    zprInit();
    zprSelectionFunc(drawSys);     /* Selection mode draw function */
    zprPickFunc(pick);              /* Pick event client callback   */

     /* Initialise OpenGL */

    GLERROR;
	
	positions = calloc(body_count*MAX_DIMS*TICK_MEMORY,sizeof(double));
	sizes = calloc(body_count,sizeof(double));
	bodyIndices = calloc(2*body_count,sizeof(int));
	
	double mass[] = {1.9891e30,1.8986e27,3.3022e23};
	double radius[] = {696342e3,69911e3,2439.7e3};
	Vector position[] = {(Vector){.e1=0,.e2=0,.e3=0},(Vector){.e1=778547200e3,.e2=0,.e3=0},(Vector){.e1=57909100e3,.e2=0,.e3=0}};
	Vector momentum[] = {(Vector){.e1=0,.e2=0,.e3=0},(Vector){.e1=0,.e2=13.07e3 * mass[1],.e3=0},(Vector){.e1=0,.e2=47.87e3 * mass[2],.e3=0}};
	
	gravsys = buildSystem(mass,radius,position,momentum,body_count);
	
	printf("System initialized \n");
	
	dt = 1000;
	t = 0;
	end_time = 3e9;//3.74336e8; //for jupiter or 7.6005e6 for mercury;
	h_init = gravsys.hamiltonian(gravsys.bodies,gravsys.n,ALL_BODIES);
	printf("Initial energy %f\n",h_init);
	
	operiod = (int)max(round(100/dt),1);
	ooperiod = TICK_MEMORY;
	ooline = oline = 0;
	for(int i = 0; i < body_count; i++){
		sizes[i] = smaller(2*radius[i]);
		bodyIndices[2*i+1] = i * MAX_DIMS * TICK_MEMORY;
	}
    /* Enter GLUT event loop */
	printf("Launching GLUT mainloop\n");
    glutMainLoop();

    return EXIT_SUCCESS;
}
