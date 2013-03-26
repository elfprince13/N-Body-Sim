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


#define LEFTH 0
#define RIGHTH 4
#define UPPERH 2
#define LOWERH 0
#define ALLQ 8

#define ZERO_VECTOR (Vector){.e1=0, .e2=0, .e3=0}
#define SET_NO_CHILDREN(node)	node->children[0][0] = NULL; \
								node->children[0][1] = NULL; \
								node->children[1][0] = NULL; \
								node->children[1][1] = NULL

#define HAS_NO_CHILDREN(node)	(node->children[0][0] == NULL && \
								 node->children[0][1] == NULL && \
								 node->children[1][0] == NULL && \
								 node->children[1][1] == NULL)

#define FREE_CHILDREN(node)	node->children[0][0] = freeTree(node->children[0][0]); \
							node->children[0][1] = freeTree(node->children[0][1]); \
							node->children[1][0] = freeTree(node->children[1][0]); \
							node->children[1][1] = freeTree(node->children[1][1])

#define MAKE_CHILDREN(node)	node->children[0][0] = initInQuad(node,LOWERH | LEFTH); \
							node->children[0][1] = initInQuad(node,UPPERH | LEFTH); \
							node->children[1][0] = initInQuad(node,LOWERH | RIGHTH); \
							node->children[1][1] = initInQuad(node,UPPERH | RIGHTH)

#define IS_DUMMY(body) (isnan(body.r))

/*
#define DOUBLEQ -1
#define HALFQ 0.5
#define SAMEQ 0
//*/

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
typedef struct{Vector lowX; Vector highX; } AABB;
typedef struct _qtNode{ AABB bounds; Body * contents; /* contents may be virtual */ struct _qtNode * parent; struct _qtNode * children [2] [2]; } QuadTree;


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
double min(double f1, double f2){ return (f1 > f2)*f2 + (f1 <= f2)*f1;	}

void copyVector(double * data, Vector * v, int * index)
{	for(int i = 0; i < MAX_DIMS; (*index)++, i++){	data[(*index)] = getVField(v,i);	}	}
double smaller(double f){ return f * 1e-10; }
void copyVectorSmaller(double * data, Vector * v, int * index)
{	Vector sml = scale(smaller(1),*v);
	copyVector(data, &sml, index);	}

AABB makeNullBB()
{
	AABB nanBB = (AABB){.lowX = (Vector){.e1 = NAN, .e2 = NAN, .e3 = NAN},
		.highX = (Vector){.e1 = NAN, .e2 = NAN, .e3 = NAN}};
	for(int i = DIMS; i < MAX_DIMS; i++) setVField(&(nanBB.lowX), i, 0),setVField(&(nanBB.highX), i, 0);
	return nanBB;
}

int intersectsAABB(AABB * bb, Vector * v)
{

	// Not necessarily centered, so much add and substract;
	int ret = 1;
	for(int i = 0; ret && i < DIMS; i++ )
		ret	= ret	&&	getVField(v,i) < getVField(&(bb->highX), i)
					&&	getVField(v,i) > getVField(&(bb->lowX), i);
	return	ret;
}

Vector dim(AABB * bb){
	Vector dim = ZERO_VECTOR;
	for(int i = 0; i < DIMS; i++ ) setVField(&dim, i, getVField(&(bb->highX), i) - getVField(&(bb->lowX), i));
	return dim;
}

Vector mid(AABB * bb)
{
	double hx, lx;
	Vector mid = ZERO_VECTOR;
	for(int i = 0; i < DIMS; i++){
		lx = getVField(&(bb->lowX), i), hx = getVField(&(bb->highX), i);
		setVField(&mid, i, lx + (hx - lx)/2);
	}
	return mid;
}

void redimAABB(AABB * bb, AABB * fixed, Vector * mag)
{
	AABB tmp = *fixed;
	double nmn,nmx,omn,omx;
	for(int i = 0; i < DIMS; i++){
		nmn = getVField(&(bb->lowX), 0);
		nmx = getVField(&(bb->highX), 0);
		omn = getVField(&(fixed->lowX), 0);
		omx = getVField(&(fixed->highX), 0);
		if(nmn == omn){
			setVField(&(tmp.lowX), i, omn);
			setVField(&(tmp.highX), i, omn+getVField(mag, i));
		} else if(nmx == omx){
				setVField(&(tmp.highX), i, omx);
				setVField(&(tmp.lowX), i, omx-getVField(mag, i));
		} else{
			printf("Error in redimAABB(): one of %G or %G should be 0\n",fabs(nmn-omn),fabs(nmx-omx));
			exit(EXIT_FAILURE);
		}
	}
	*bb = tmp;
}

void updateAABB(AABB * bb, Vector * v)
{
	double x,l,h;
	for(int i = 0; i < DIMS; i++)
	{
		x = getVField(v, i);
		l = getVField(&(bb->lowX), i);
		h = getVField(&(bb->highX), i);
		if(x < l || isnan(l)) setVField(&(bb->lowX), i, x - sqEP);
		if(x > l || isnan(h)) setVField(&(bb->highX), i, x + sqEP);
	}
}

Vector barycenter(Body * b1, Body * b2)
{
	return scale((1/(b1->m + b2->m)),sum(scale(b1->m, b1->x),scale(b2->m, b2->x)));
}

Vector corner(Vector *v1, Vector *v2, int quad)
{
	Vector corner = ZERO_VECTOR;
	for(int i = 0; i < DIMS; i++){
		setVField(&corner, i, ((quad >> (MAX_DIMS - i - 1)) & 0x1) ?	max(getVField(v1, i),getVField(v2, i)) :
																		min(getVField(v1, i),getVField(v2, i)));
	}
	return corner;
}

void updateMasses(QuadTree * node, Body * b)
{
	
	if(node->contents != NULL)
	{
		Body *ob = node->contents;
		Vector nbc = barycenter(ob,b);
		double nm = ob->m + b->m;
		
		// Don't *ever* do this to a real body
		// Only virtual bodies for interior nodes
		ob->x = nbc;
		ob->m = nm;
	} else {
		node->contents = b;
	}
}

Body* getDummyBody()
{
	Body * b = calloc(1,sizeof(Body));
	b->m = 0;
	b->r = NAN;
	b->x = ZERO_VECTOR;
	b->p = ZERO_VECTOR;
	return b;
}

void returnDummyBody(QuadTree * qt)
{
	// If this happens a lot, make a dummy-body pool
	if(qt->contents != NULL && IS_DUMMY( (*(qt->contents)) )){
		free(qt->contents);
		qt->contents = NULL;
	}
}

QuadTree * freeTree(QuadTree * node)
{
	returnDummyBody(node);
	FREE_CHILDREN(node);
	free(node);
	return NULL;
}

QuadTree * initEmpty()
{
	QuadTree * ret = calloc(1, sizeof(QuadTree));
	SET_NO_CHILDREN(ret);
	ret->parent = NULL;
	ret->bounds = makeNullBB();
	ret->contents = NULL;
	return ret;
	
}

QuadTree * initInQuad(QuadTree * node, int quad)
{
	QuadTree * ret = initEmpty();
	ret->parent = node;
	Vector m = mid(&(node->bounds));
	Vector c = ZERO_VECTOR;
	Vector lx = ZERO_VECTOR;
	Vector hx = ZERO_VECTOR;
	for(int i = 0; i < DIMS; i++){
		setVField(&c, i, getVField(((quad >> (MAX_DIMS - i - 1)) & 0x1) ? &(node->bounds.highX) : &(node->bounds.lowX),i));
		setVField(&lx, i, min(getVField(&m, i),getVField(&c, i)));
		setVField(&hx, i, max(getVField(&m, i),getVField(&c, i)));
	}
	return ret;
	
	
}

void initChildren(QuadTree * node)
{
	MAKE_CHILDREN(node);
	
}

void splitNodeB(QuadTree * node, Body * b)
{
	if(HAS_NO_CHILDREN(node))
	{
		if(node->contents == NULL){
			node->contents = b;
		} else{
			initChildren(node);
			
		}
	}
	if(node->contents != b)
	{
		
	}
}

void splitNodeQ(QuadTree * node, QuadTree * descend)
{
	
}

// Build a region-bounding quad-tree bottom up
// More complex branching
// But saves a lead-in time of N building initial AABB
QuadTree * insertBody(QuadTree * node, Body * b)
{
	if(node == NULL)
	{
		node = calloc(1, sizeof(QuadTree));
		node->parent = NULL;
		SET_NO_CHILDREN(node);
		
		node->bounds = makeNullBB();
		
		node->contents = b;
		updateAABB(&(node->bounds), &(b->x));
	} else if(node->parent == NULL && HAS_NO_CHILDREN(node)) {
		// Nobody here but us chickens
		updateAABB(&(node->bounds), &(b->x));
		
		// Split Node
		splitNodeB(node, b);
		
	} else if(intersectsAABB(&(node->bounds), &(b->x))){
		// Add to a child or split
		if(HAS_NO_CHILDREN(node))
		{
			splitNodeB(node, b);
		} else{
			// Safe to discard return value,
			// since safely inside!
			Vector nmid = mid(&(node->bounds));
			insertBody(node->children	[(b->x.e1 > nmid.e1)]
										[(b->x.e2 > nmid.e2)], b);
		}
	} else {
		// reparent the tree
		Vector cdim = dim(&(node->bounds));
		Vector  ndim, dimr;
		double nsubd;
		QuadTree * nrent = calloc(1, sizeof(QuadTree));
		SET_NO_CHILDREN(nrent);
		nrent->bounds = node->bounds;
		updateAABB(&(nrent->bounds), &(b->x));
		
		ndim =dim(&(nrent->bounds));
		nsubd = max(ceil(log2(ndim.e1 / cdim.e1)),ceil(log2(ndim.e2 / cdim.e2)));
		dimr = scale(pow(2,nsubd),cdim);
		
		redimAABB(&(nrent->bounds), &(node->bounds), &dimr);
		splitNodeQ(nrent, node);
		node = nrent;
		
	}
	
	return node;
}


Vector gradient(State s,int n,int kind)
{
	Vector grad = ZERO_VECTOR;
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

void drawTree(QuadTree * node,int depth)
{
	static double verts[MAX_DIMS * 4];
	static unsigned char indices[8] = {0,1,1,2,2,3,3,0};
	static float colors[3] = {0.9,0.3,0.6};
	int i = 0;
	
	if(node == NULL) return;
	
	Vector corners[2][2] = {{	corner(&(node->bounds.lowX),&(node->bounds.highX),LOWERH | LEFTH ),
								corner(&(node->bounds.lowX),&(node->bounds.highX),UPPERH | LEFTH )},
							{	corner(&(node->bounds.lowX),&(node->bounds.highX),LOWERH | RIGHTH ),
								corner(&(node->bounds.lowX),&(node->bounds.highX),UPPERH | RIGHTH )}};
	copyVector(verts, &corners[0][0], &i);
	copyVector(verts, &corners[0][1], &i);
	copyVector(verts, &corners[1][0], &i);
	copyVector(verts, &corners[1][1], &i);
	
	glDisable(GL_CULL_FACE);
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(MAX_DIMS*4,GL_DOUBLE,0,verts);
	
	glColor3f(colors[depth%3],colors[(depth+1)%3],colors[(depth+2)%3]);
	glLineWidth(1);
	glDrawElements(GL_LINES, 4, GL_UNSIGNED_BYTE, indices);
	
	glDisableClientState(GL_VERTEX_ARRAY);
	
	
	glEnable(GL_CULL_FACE);
	drawTree(node->children[0][0],depth+1);
	drawTree(node->children[0][1],depth+1);
	drawTree(node->children[1][0],depth+1);
	drawTree(node->children[1][1],depth+1);
	
}

void drawTime()
{
	static char buf[128];
	sprintf(buf,"Time: %G / %G with %.03G%% H_drift",t,end_time,(h_drift / h_init) * 100);
	glRasterPos2i(-100,100);
	glColor3f(1.0f, 1.0f, 0.3f);
	for(char *b=buf;*b;b++) glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *b);
}

/* Callback functions for drawing */
static QuadTree * testnode;
void tdisplay(int done)
{
   GLERROR;

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// printf("Stepping time from %f\n",t);
	stepSys(TICKS_PER);
	//printf("Drawing system\n");
	drawSys();
	drawTree(testnode, 0);
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
	Vector position[] = {ZERO_VECTOR,(Vector){.e1=778547200e3,.e2=0,.e3=0},(Vector){.e1=57909100e3,.e2=0,.e3=0}};
	Vector momentum[] = {ZERO_VECTOR,(Vector){.e1=0,.e2=13.07e3 * mass[1],.e3=0},(Vector){.e1=0,.e2=47.87e3 * mass[2],.e3=0}};
	
	gravsys = buildSystem(mass,radius,position,momentum,body_count);
	
	printf("System initialized \n");
	
	testnode = insertBody(NULL, gravsys.bodies[0]);
	testnode = insertBody(testnode, gravsys.bodies[1]);
	
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
