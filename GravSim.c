#include <stdio.h>
#include <stdlib.h>

#include "GravSim.h"
//#include "NaiveIntegrator.h"
#include "RungeKutta.h"
#include "QuadTree.h"
#include "QTNodeList.h"
#include "Constants.h"

static double _matrix[16];
GLdouble * getMatrix()
{
	glGetDoublev(GL_MODELVIEW_MATRIX,_matrix);
	return _matrix;
}

double determinant4x4(const GLdouble *m)
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
#undef MAT
#undef m11
#undef m12
#undef m13
#undef m14
#undef m21
#undef m22
#undef m23
#undef m24
#undef m31
#undef m32
#undef m33
#undef m34
#undef m41
#undef m42
#undef m43
#undef m44
}


//--------------------------------------------------
// System parameters
//--------------------------------------------------
static State gravsys;
static NGradient gradient;
/*static*/ int body_count = 11;

static double theta = 0.85;
static QuadTree * primaryBHTree=NULL;
int useBarnes = 0;

/*static*/ double dt = 10;//.1;
static double t = 0;
static double end_time = /*100000000;*/ 3.74336e8; //for jupiter or 7.6005e6 for mercury;
static double softening = 5000000000;

static int operiod;
static int  oline;
static int ooperiod;
static int ooline;

static double h_init, h_drift;

static double * positions;
static double * sizes;
static int * bodyIndices;

static int follow = -1;
static double fzoom = 1;


double tc_start = 0.0;
double tc_end = 0.0;
long frames = 0;
char titlestring[200];
double fps;
static int should_auto_exit;


void showFPS() {
	
    /*double t;
    
    // Get current time
    t = .001*glutGet(GLUT_ELAPSED_TIME);  // Get number of seconds since glfwInit()
						// Calculate and display FPS (frames per second) in title bar.
    if( (t-t0) > 0.5 || frames == 0 )
    {*/
	fps = TICKS_PER / (tc_end - tc_start);//(double)frames / (t-t0);
	//   t0 = t;
    //    frames = 0;
    //}
    frames ++;
	sprintf(titlestring, "(%.1f FPS)", fps);
	glRasterPos2i(100,-100);
	glColor3f(1.0f, 1.0f, 0.3f);
	for(char *b=titlestring;*b;b++) glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *b);
	//printf("%f\n",fps);
}

static RKMethod rkMethod = (RKMethod){
	.tableau = NULL,
	.init_tableau = NULL,
	.integrate = NULL,
	.stages = 0
};

QuadTree * barnesHutIterator(State * s,QuadTree *root,int which){
	static State * internalState = NULL;
	static int n = -1;
	static Body * curBody = NULL;
	static QTList * todo = NULL;
	static QTList * outl = NULL;
	static QuadTree * internalRoot = NULL;
	
	if(which != n){
		freeQTList(todo);
		freeQTList(outl);
		curBody = NULL;
		internalRoot = root;
		internalState = s;
		n = which;
		if(n != -1){
			curBody = internalState->bodies[n];
			todo = (QTList*)malloc(sizeof(QTList));
			todo->head = NULL;
			todo->tail = NULL;
			outl = (QTList*)malloc(sizeof(QTList));
			outl->head = NULL;
			outl->tail = NULL;
			append(todo, internalRoot);
		} else{
			todo = NULL;
			outl = NULL;
			return NULL;
		}
	}
	
	QuadTree * cur;
	while(!isEmptyQTList(todo)){
		cur = pophead(todo);
		if(cur->contents == NULL)
			continue;
		
		Vector boxSize = diff(cur->bounds.highX,cur->bounds.lowX);
		double l = max(boxSize.e1,boxSize.e2);
		double d = mag(diff(curBody->x,cur->contents->x));
		if(/*intersectsAABB(&(cur->bounds), &(curBody->x)) ||*/ d <= l/theta){
			if(!HAS_NO_CHILDREN(cur)){
				append(todo,cur->children[0][0]);
				append(todo,cur->children[0][1]);
				append(todo,cur->children[1][0]);
				append(todo,cur->children[1][1]);
			} else if(cur->contents != NULL && !intersectsAABB(&(cur->bounds), &(curBody->x))){
				append(outl,cur);
				break;
			}
		} else{
			append(outl,cur);
			break;
		}
	}
	
	return pophead(outl);
}



Vector numericalGradient(State * s,int n,int kind)
{
	Vector grad = ZERO_VECTOR;
	double H_l,H_h,h,e;
	Body ** bodies = s->bodies;
	Body * b = bodies[n];
	Body midpoint = *b;
	Body m_l, m_h;
	double h0 = s->hamiltonian(s,n);
	for(int i = kind; i < kind + DIMS; i++){
		h = /**/fabs(sqEP * (sqEP + getBField(b,i)));/**//*sqEP * sqrt(sqrt(sqEP+fabs(h0))*fabs(sqEP + getBField(b,i)))*/;
		m_l = midpoint;
		m_h = midpoint;
		difBField(&m_l,i,-h);
		difBField(&m_h,i,h);
		bodies[n] = &m_l;
		
		H_l = s->hamiltonian(s,n);
		bodies[n] = &m_h;
		H_h = s->hamiltonian(s,n);
		
		setVField(&grad,i%MAX_DIMS,(H_h - H_l) / (2*h));
		bodies[n] = b;
	}
	return grad;
}

// The selector allows us to determine whether we want a single body or the whole system
// This should speed up differentiation by an order of N, since we may neglect terms
// which will not vary.
double newtonianGravitation(State * s, int which){
	double H = 0;
	double H_i = 0;
	
	Body ** bodies = s->bodies;
	int n = s->n;
	
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

// The selector allows us to determine whether we want a single body or the whole system
// This should speed up differentiation by an order of N, since we may neglect terms
// which will not vary.
double softNewtonianGravitation(State * s, int which){
	double H = 0;
	double H_i = 0;
	
	Body ** bodies = s->bodies;
	int n = s->n;
	
	int istart = (which < 0) ? 0 : which;
	int iend = (which < 0) ? n : which + 1;
	int jend, i, j;
	for(i = istart; i < iend; i++){
		H_i = 0;
		jend = (which < 0) ? i : n;
		for(j  = !i ? 1 : 0; j < jend; j = (j+1 == i) ? j + 2 : j + 1){
			H_i += bodies[i]->m * bodies[j]->m / (softening+ mag(diff(bodies[i]->x,bodies[j]->x)));
		}
		H_i = dot(bodies[i]->p,bodies[i]->p) / (2 * bodies[i]->m) - G*(H_i);
		H += H_i;
	}
	return H;
}

// The selector allows us to determine whether we want a single body or the whole system
// This should speed up differentiation by an order of N, since we may neglect terms
// which will not vary.
double barnesHutGravitation(State * s, int which){
	double H = 0;
	double H_i = 0;
	
	Body ** bodies = s->bodies;
	QuadTree * bhN = NULL;
	int n = s->n;
	double wfactor = (which < 0 ) ? 2 : 1;
	
	int istart = (which < 0) ? 0 : which;
	int iend = (which < 0) ? n : which + 1;
	int i;
	for(i = istart; i < iend; i++){
		H_i = 0;
		while((bhN = barnesHutIterator(s, primaryBHTree, i)) != NULL){
			H_i += bodies[i]->m * bhN->contents->m / mag(diff(bodies[i]->x,bhN->contents->x));
		}
		H_i = dot(bodies[i]->p,bodies[i]->p) / (2 * bodies[i]->m) - G*(H_i) / wfactor;
		H += H_i;
	}
	barnesHutIterator(NULL,NULL,-1);
	return H;
}

// The selector allows us to determine whether we want a single body or the whole system
// This should speed up differentiation by an order of N, since we may neglect terms
// which will not vary.
double softBarnesHutGravitation(State * s, int which){
	double H = 0;
	double H_i = 0;
	
	Body ** bodies = s->bodies;
	QuadTree * bhN = NULL;
	int n = s->n;
	double wfactor = (which < 0 ) ? 2 : 1;
	
	int istart = (which < 0) ? 0 : which;
	int iend = (which < 0) ? n : which + 1;
	int i;
	for(i = istart; i < iend; i++){
		H_i = 0;
		while((bhN = barnesHutIterator(s, primaryBHTree, i)) != NULL){
			H_i += bodies[i]->m * bhN->contents->m / (softening + mag(diff(bodies[i]->x,bhN->contents->x)));
		}
		H_i = dot(bodies[i]->p,bodies[i]->p) / (2 * bodies[i]->m) - G*(H_i) / wfactor;
		H += H_i;
	}
	barnesHutIterator(NULL,NULL,-1);
	return H;
}

Vector newtonianGravitationGradient(State * s, int which, int kind){
	Body ** bodies = s->bodies;
	int n = s->n;
	int i = which;
	int j;
	Vector ret,gvec,lvec,dvec;
	double inmag,gm;
	gm = bodies[i]->m;
	switch(kind){
		case X_START:
			ret = ZERO_VECTOR;
			gvec = bodies[i]->x;
			for (j  = !i ? 1 : 0; j < n; j = (j+1 == i) ? j + 2 : j + 1)
			{
				lvec = bodies[j]->x;
				dvec = diff(lvec,gvec);
				inmag = 1/mag(dvec);
				ret = sum(ret,scale(gm*(bodies[j]->m)*inmag*inmag*inmag,dvec));
			}
			ret = scale(-G,ret);
			break;
		case P_START:
			gvec = bodies[i]->p;
			ret = scale(1/(gm),gvec);
			break;
		default:
			printf("Error in newtonianGravitationGradient():  %d is not a valid specifier for kind\n",kind);
	}
	return ret;
}

Vector softNewtonianGravitationGradient(State * s, int which, int kind){
	Body ** bodies = s->bodies;
	int n = s->n;
	int i = which;
	int j;
	Vector ret,gvec,lvec,dvec;
	double inmag,gm;
	gm = bodies[i]->m;
	switch(kind){
		case X_START:
			ret = ZERO_VECTOR;
			gvec = bodies[i]->x;
			for (j  = !i ? 1 : 0; j < n; j = (j+1 == i) ? j + 2 : j + 1)
			{
				lvec = bodies[j]->x;
				dvec = diff(lvec,gvec);
				inmag = 1/(softening+mag(dvec));
				ret = sum(ret,scale(gm*(bodies[j]->m)*inmag*inmag*inmag,dvec));
			}
			ret = scale(-G,ret);
			break;
		case P_START:
			gvec = bodies[i]->p;
			ret = scale(1/(gm),gvec);
			break;
		default:
			printf("Error in newtonianGravitationGradient():  %d is not a valid specifier for kind\n",kind);
	}
	return ret;
}

Vector barnesHutGravitationGradient(State * s, int which, int kind){
	Body ** bodies = s->bodies;
	int n = s->n;
	int i = which;
	QuadTree * bhN;
	Vector ret,gvec,lvec,dvec;
	double inmag,gm;
	gm = bodies[i]->m;
	switch(kind){
		case X_START:
			ret = ZERO_VECTOR;
			gvec = bodies[i]->x;
			while((bhN = barnesHutIterator(s, primaryBHTree, i)) != NULL){
				lvec = bhN->contents->x;
				dvec = diff(lvec,gvec);
				inmag = 1/mag(dvec);
				ret = sum(ret,scale(gm*(bhN->contents->m)*inmag*inmag*inmag,dvec));
			}
			barnesHutIterator(NULL,NULL,-1);
			ret = scale(-G,ret);
			break;
		case P_START:
			gvec = bodies[i]->p;
			ret = scale(1/(gm),gvec);
			break;
		default:
			printf("Error in newtonianGravitationGradient():  %d is not a valid specifier for kind\n",kind);
	}
	return ret;
}

Vector softBarnesHutGravitationGradient(State * s, int which, int kind){
	Body ** bodies = s->bodies;
	int n = s->n;
	int i = which;
	QuadTree * bhN;
	Vector ret,gvec,lvec,dvec;
	double inmag,gm;
	gm = bodies[i]->m;
	switch(kind){
		case X_START:
			ret = ZERO_VECTOR;
			gvec = bodies[i]->x;
			while((bhN = barnesHutIterator(s, primaryBHTree, i)) != NULL){
				lvec = bhN->contents->x;
				dvec = diff(lvec,gvec);
				inmag = 1/(softening+mag(dvec));
				ret = sum(ret,scale(gm*(bhN->contents->m)*inmag*inmag*inmag,dvec));
			}
			barnesHutIterator(NULL,NULL,-1);
			ret = scale(-G,ret);
			break;
		case P_START:
			gvec = bodies[i]->p;
			ret = scale(1/(gm),gvec);
			break;
		default:
			printf("Error in newtonianGravitationGradient():  %d is not a valid specifier for kind\n",kind);
	}
	return ret;
}

Body bDeriv(State * s, int which, Body * testpos, int kind)
{
	Body ret = (Body){.x = ZERO_VECTOR, .p = ZERO_VECTOR, .m = testpos->m, .r = testpos-> r};
	Body * opos = s->bodies[which];
	s->bodies[which] = testpos;
	
	ret.x = (kind == X_START || kind == X_AND_P) ? gradient(s,which,P_START) : ZERO_VECTOR;
	ret.p = (kind == P_START || kind == X_AND_P) ? scale(-1,gradient(s,which,X_START)) : ZERO_VECTOR;

	s->bodies[which] = opos;
	return ret;
}

State buildSystem(NHamiltonian nh, NGradient ng, double m[], double r[], Vector x[], Vector p[], int n){
	Body ** bodies = (Body**)calloc(sizeof(Body*),n);
	for(int i = 0; i < n; i++)	*(bodies[i] = (Body*)calloc(sizeof(Body),1)) = (Body){.x=x[i],.p=p[i],.m=m[i],.r=r[i]};
	return (State){.hamiltonian = nh,.analyticalGradient = ng,.n = n,.bodies = bodies};
}

int teardownSystem(State s){
	for(int i = 0; i < s.n; i++) free(s.bodies[i]);
	free(s.bodies);
	s.bodies = NULL;
	return 0;
}

void setIntegrator(TableauGen);

void stepSys(int n){
	double h_tot;
	//double h_tot2;
	Body deltaBodies[body_count];
	//Vector xdot[body_count],pdot[body_count];
	//char buf[64][64];
	for(;t < end_time && n > 0; t += dt, oline = (oline + 1) % operiod, ooline = (ooline + !oline) % ooperiod, n--){
		if(useBarnes){
			primaryBHTree = initEmptyQuad();
			for(int i = 0; i < body_count; updateAABB(&(primaryBHTree->bounds), &(gravsys.bodies[i++]->x)));
			for(int i = 0; i < body_count; insertBody(primaryBHTree, gravsys.bodies[i++]));
		}
		
		h_tot = gravsys.hamiltonian(&gravsys,ALL_BODIES);
		//h_tot2 = newtonianGravitation(&gravsys, ALL_BODIES);
		h_drift = h_init - h_tot;
		if(!oline){
			/*
			 printf("%G %G %G %s %s %s %s %s %s\n",t,h_tot,h_drift,
						stringify(gravsys.bodies[0]->x,buf[0]),stringify(gravsys.bodies[0]->p,buf[1]),
						stringify(gravsys.bodies[1]->x,buf[2]),stringify(gravsys.bodies[1]->p,buf[3]),
						stringify(gravsys.bodies[2]->x,buf[4]),stringify(gravsys.bodies[2]->p,buf[5]));
			 //*/
			
			
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
			
			/*setIntegrator(init_PRK6);
			Body b1 = rkMethod.integrate(&gravsys,i,bDeriv,&rkMethod);
			setIntegrator(init_RK4);
			Body b2 = rkMethod.integrate(&gravsys,i,bDeriv,&rkMethod);
			setIntegrator(init_ForwardEuler);
			Body b3 = rkMethod.integrate(&gravsys,i,bDeriv,&rkMethod); */
			//printf("%X %X %X",&b1, &b2,&b3);
			/*
			Vector xdota,xdotn,pdota,pdotn;
			xdotn = numericalGradient(&gravsys, i, X_START);
			xdota = newtonianGravitationGradient(&gravsys, i, X_START);
			pdotn = numericalGradient(&gravsys, i, P_START);
			pdota = newtonianGravitationGradient(&gravsys, i, P_START); */
			deltaBodies[i] = rkMethod.integrate(&gravsys,i,bDeriv,&rkMethod);

		}
		if(useBarnes){
			freeTree(primaryBHTree);
			primaryBHTree=NULL;
		}
		for(int i = 0; i < body_count; i++){
			for(int j = 0; j < DIMS; j++){
				difBField(gravsys.bodies[i],X_START+j,getBField(&deltaBodies[i],X_START+j));
				difBField(gravsys.bodies[i],P_START+j,getBField(&deltaBodies[i],P_START+j));
			}
		}
	}
}

#define GLERROR                                                    \
    {                                                              \
        GLenum code = glGetError();                                \
        while (code!=GL_NO_ERROR)                                  \
        {                                                          \
            printf("%s (%d)\n",(char *) gluErrorString(code), code);          \
                code = glGetError();                               \
        }                                                          \
    }

void drawSys()
{
	static float colors[3] = {0.9,0.3,0.6};
	double scale = pow(determinant4x4(getMatrix()),1/3.);
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
		glPointSize(5);//sizes[i]*scale);
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
	static unsigned char indices[/*8*/5] = {0,1,2,3,0};//{0,1,1,2,2,3,3,0};
	static float colors[10] = {0.9,0.,0.6,0.3,0.3,0.6,0.,0.9,0.3,0.6};
	int i = 0;
	
	
	if(node == NULL) return;
	
	Vector corners[2][2] = {{	corner(&(node->bounds.lowX),&(node->bounds.highX),LOWERH | LEFTH ),
								corner(&(node->bounds.lowX),&(node->bounds.highX),UPPERH | LEFTH )},
							{	corner(&(node->bounds.lowX),&(node->bounds.highX),LOWERH | RIGHTH ),
								corner(&(node->bounds.lowX),&(node->bounds.highX),UPPERH | RIGHTH )}};
	copyVectorSmaller(verts, &corners[0][0], &i);
	copyVectorSmaller(verts, &corners[0][1], &i);
	copyVectorSmaller(verts, &corners[1][1], &i);
	copyVectorSmaller(verts, &corners[1][0], &i);
	
	glDisable(GL_CULL_FACE);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3,GL_DOUBLE,0,verts);
	glColor3f(colors[depth%10],colors[(depth+1)%10],colors[(depth+2)%10]);
	glLineWidth(3);
	glDrawElements(GL_LINE_STRIP, 5, GL_UNSIGNED_BYTE, indices);
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
void tdisplay(int done)
{
   GLERROR;

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// printf("Stepping time from %f\n",t);
	tc_start = .001*glutGet(GLUT_ELAPSED_TIME);
	stepSys(TICKS_PER);
	tc_end = .001*glutGet(GLUT_ELAPSED_TIME);

	if(follow > -1){
		//getMatrix();
		glLoadIdentity();
		glScalef(.008,.008,.008);
		Vector x = gravsys.bodies[follow]->x;
		glTranslated(-smaller(x.e1)*fzoom,-smaller(x.e2)*fzoom,-smaller(x.e3)*fzoom);
		glScalef(fzoom,fzoom,fzoom);
		
	}
	//printf("Drawing system\n");
	drawSys();
	
	//*
	primaryBHTree = initEmptyQuad();
	for(int i = 0; i < body_count; updateAABB(&(primaryBHTree->bounds), &(gravsys.bodies[i++]->x)));
	for(int i = 0; i < body_count; insertBody(primaryBHTree, gravsys.bodies[i++]));
	drawTree(primaryBHTree, 0);
	freeTree(primaryBHTree);
	primaryBHTree=NULL;
	//*/
	drawTime();
	showFPS();
	//if(!done && t<=end_time && fps != INFINITY) printf("%G %G %G %G %f\n",t,h_init,h_init+h_drift,h_drift,fps);

	//printf("next frame\n\n");
   glutSwapBuffers();

   GLERROR;
	if(!done) glutTimerFunc(17,tdisplay,t >= end_time);
	else if(done && t>=end_time && should_auto_exit) exit(0);
}
void display(){ tdisplay(1); /* Don't schedule extra timers */ }

/* Callback function for pick-event handling from ZPR */

void pick(GLint name)
{
	//printf("Pick: %d\n",name);
	//fflush(stdout);
}


void setIntegrator(TableauGen name){
	
	double **a, **A, *b,*B,*c,*C,*d;
	
	if(rkMethod.init_tableau == init_ForwardEuler || rkMethod.init_tableau == init_RK4){
		((StandardTableau*)rkMethod.tableau)->tfree(rkMethod.tableau);
	}
	if(rkMethod.init_tableau == init_GLRK4){
		((NonSingularImplicitTableau*)rkMethod.tableau)->tfree(rkMethod.tableau);
	}
	if(rkMethod.init_tableau == init_PRK6Ruth || rkMethod.init_tableau == init_PRK6SuzukiTrotter || rkMethod.init_tableau == init_PRK6BlanesMoan){
		((PartitionedTableau*)rkMethod.tableau)->tfree(rkMethod.tableau);
	}
	
	int stages;
	if(name == init_ForwardEuler || name == init_RK4){
		stages = (name == init_ForwardEuler) ? 1 : 4;
		StandardTableau * st = malloc(sizeof(StandardTableau));
		
		a = makeDynamic2DArray(stages, stages);
		b = malloc(stages*sizeof(double));
		c = malloc(stages*sizeof(double));
		
		*st = (StandardTableau){.a = a, .b = b, .c = c, .tfree = free_StandardRK};
		
		rkMethod = (RKMethod){
			.init_tableau = name,
			.tableau = (void*)(st),
			.integrate = integrate_RK,
			.stages = stages,
		};
	} else if(name == init_PRK6Ruth || name == init_PRK6SuzukiTrotter || name == init_PRK6BlanesMoan){
		stages = 6;
		PartitionedTableau * st = malloc(sizeof(PartitionedTableau));
		
		c = malloc(stages*sizeof(double));//makeDynamic2DArray(stages, stages);
		b = malloc(stages*sizeof(double));
		/*c = malloc(stages*sizeof(double));
		A = makeDynamic2DArray(stages, stages);
		B = malloc(stages*sizeof(double));
		C = malloc(stages*sizeof(double));*/
		*st = (PartitionedTableau){.a = /*a*/c, .b = b,/* .c = c,.A = A, .B = B, .C = C, */.tfree = free_PartitionedRK};
		rkMethod = (RKMethod){
			.init_tableau = name,
			.tableau = (void*)(st),
			.integrate = integrate_PRK,
			.stages = stages
		};
		
	} else {
		rkMethod = (RKMethod){
			.init_tableau = NULL,
			.tableau = NULL,
			.integrate = NULL,
			.stages = 0,
		};
	}
	
	if(rkMethod.init_tableau != NULL){
		rkMethod.init_tableau(rkMethod.tableau);
	}
}

void KeyHandler(unsigned char key, int x, int y){
	switch(key){
		case ' ':
			follow = (follow+2)%(body_count+1) - 1;
			if(follow < 0) fzoom = 1;
			printf("following %d\n",follow);
			break;
		case '-':
			fzoom /= 1.4;
			break;
		case '=':
			fzoom *= 1.4;
		default:
			0;//printf("%c was pressed\n",key);
	}
}

/* Entry point */

int main(int argc, char *argv[])
{
	TableauGen integrator_select[4] = {init_ForwardEuler, init_RK4, init_PRK6Ruth, init_PRK6SuzukiTrotter};
	NHamiltonian hamilton_select[4] = {newtonianGravitation, barnesHutGravitation,
										softNewtonianGravitation,softBarnesHutGravitation};
	NGradient gradient_select[4] = {newtonianGravitationGradient, barnesHutGravitationGradient,
									softNewtonianGravitationGradient,softBarnesHutGravitationGradient};
	int body_setup;
	double max_time_select, dt_select;
	int auto_exit_select;
	
	int iselector,hselector;
	
	if(argc != 7 && argc != 1){
		printf("Usage: %s body_setup dt max_time auto_exit integrator hamiltonian_impl\n",argv[0]);
		exit(1);
	} else if (argc == 7){
		body_setup = atoi(argv[1]);
		dt_select = atof(argv[2]);
		max_time_select = atof(argv[3]);
		auto_exit_select = atoi(argv[4]);
		iselector = atoi(argv[5]);
		hselector = atoi(argv[6]);
		//printf("Selecting %d\n",hselector);
	} else{
		body_setup  = 0;
		dt_select = 1000;
		max_time_select = 3e9;
		auto_exit_select = 0;
		iselector = 0;
		hselector = 0;
	}
	
	useBarnes = hselector == 1 || hselector == 3;
	
	setIntegrator(integrator_select[iselector]);
	should_auto_exit = auto_exit_select;
	
    /* Initialise GLUT and create a window */

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800,600);
    glutCreateWindow("GravSim");

    /* Configure GLUT callback functions */

    glutDisplayFunc(display);
	glutTimerFunc(17,tdisplay,t >= end_time);
	glutKeyboardFunc(KeyHandler);

    glScalef(.008,.008,.008);

    /* Configure ZPR module */

    zprInit();
    zprSelectionFunc(drawSys);     /* Selection mode draw function */
    zprPickFunc(pick);              /* Pick event client callback   */

     /* Initialise OpenGL */

    GLERROR;
	body_count = body_setup ? body_setup : 11;
	
	positions = calloc(body_count*MAX_DIMS*TICK_MEMORY,sizeof(double));
	sizes = calloc(body_count,sizeof(double));
	bodyIndices = calloc(2*body_count,sizeof(int));
	
	double * mass, * radius;
	Vector * position, * momentum;
	
	if(!body_setup){
	// Sun, Jupiter,
	// Io,Callisto,Europa,Ganymede,
	// Mercury, Venus, Earth, Mars
	
	double ssmass[] = {1.9891e30,1.8986e27,
		8.9319e22,1.075938e23,4.7998e22,1.4819e23, //Jovian moons
		3.3022e23, 4.8685e24,5.9736e24,7.3477e22, 6.4185e23}; // Inner planets and Earth-moon
	double ssradius[] = {696342e3,69911e3,
		1821.3e3,2410.3e3,1560.8e3,2634.1e3, // Jovian moons
		2439.7e3,6051.8e3,6371e3,1737.1e3,3386e3}; /// Inner planets and Earth-moon
	Vector ssposition[] = {ZERO_VECTOR,(Vector){.e1=778547200e3,.e2=0,.e3=0},
		// Jovian Moons
		(Vector){.e1=778547200e3+421700e3,.e2=0,.e3=0},
		(Vector){.e1=778547200e3+1882700e3,.e2=0,.e3=0},
		(Vector){.e1=778547200e3+670900e3,.e2=0,.e3=0},
		(Vector){.e1=778547200e3+1070400e3,.e2=0,.e3=0},
		// Inner Planets and Earth-moon
		(Vector){.e1=57909100e3,.e2=0,.e3=0},
		(Vector){.e1=108208000e3,.e2=0,.e3=0},
		(Vector){.e1=149598261e3,.e2=0,.e3=0},
		(Vector){.e1=149598261e3+384399e3,.e2=0,.e3=0},
		
		
		(Vector){.e1=227939100e3,.e2=0,.e3=0}};
	Vector ssmomentum[] = {ZERO_VECTOR,(Vector){.e1=0,.e2=13.07e3 * ssmass[1],.e3=0},
		// Jovian Moons
		(Vector){.e1=0,.e2=(13.07e3+17.334e3) * ssmass[2],.e3=0},
		(Vector){.e1=0,.e2=(13.07e3+8.204e3) * ssmass[3],.e3=0},
		(Vector){.e1=0,.e2=(13.07e3+13.740e3) * ssmass[4],.e3=0},
		(Vector){.e1=0,.e2=(13.07e3+10.880e3) * ssmass[5],.e3=0},
		// Inner Planets and Earth-moon
		(Vector){.e1=0,.e2=47.87e3 * ssmass[6],.e3=0},
		(Vector){.e1=0,.e2=35.02e3 * ssmass[7],.e3=0},
		(Vector){.e1=0,.e2=29.78e3 * ssmass[8],.e3=0},
		(Vector){.e1=0,.e2=(29.78e3+1.022e3) * ssmass[9],.e3=0},
	
		(Vector){.e1=0,.e2=24.077e3 * ssmass[10],.e3=0}};
		
		mass = ssmass;
		radius = ssradius;
		position = ssposition;
		momentum = ssmomentum;
	} else{
		srand(body_count);
		// This will leak. Oh well \_0_/
		mass = calloc(body_count,sizeof(double));
		radius = calloc(body_count,sizeof(double));
		position = calloc(body_count,sizeof(Vector));
		momentum = calloc(body_count,sizeof(Vector));
		for(int i = 0; i < body_count; i++){
			mass[i] = (rand()%10+1)*pow(10,rand()%7 + 23);
			radius[i] = pow(mass[i],0.33333)/(rand()%5);
			position[i] = (Vector){.e1=(rand()%2000-1000)*pow(10,rand()%2 + 8),.e2=(rand()%2000-1000)*pow(10,rand()%2 + 8),.e3=0};
			momentum[i] = (Vector){.e1=(rand()%60-30)*mass[i],.e2=(rand()%60-30)*mass[i],.e3=0};
		}
	}
	//printf("Now...\n");
	gravsys = buildSystem(hamilton_select[hselector],gradient_select[hselector], mass,radius,position,momentum,body_count);
	//printf("Selected %d\n",hselector);
	gradient = (gravsys.analyticalGradient == NULL) ? &numericalGradient : gravsys.analyticalGradient ;
	
	//printf("System initialized \n");
	
	
	dt = dt_select;
	t = 0;
	end_time = max_time_select;//3.74336e8; //for jupiter or 7.6005e6 for mercury;
	if(useBarnes){
		primaryBHTree = initEmptyQuad();
		for(int i = 0; i < body_count; updateAABB(&(primaryBHTree->bounds), &(gravsys.bodies[i++]->x)));
		for(int i = 0; i < body_count; insertBody(primaryBHTree, gravsys.bodies[i++]));
	}
	h_init = gravsys.hamiltonian(&gravsys,ALL_BODIES);
	if(useBarnes) freeTree(primaryBHTree);
	primaryBHTree=NULL;
	printf("%G %G %G %G %f\n",t,h_init,h_init+h_drift,h_drift,fps);
	
	operiod = (int)max(round(100/dt),1);
	ooperiod = TICK_MEMORY;
	ooline = oline = 0;
	for(int i = 0; i < body_count; i++){
		sizes[i] = smaller(2*radius[i]);
		bodyIndices[2*i+1] = i * MAX_DIMS * TICK_MEMORY;
	}
    /* Enter GLUT event loop */
	//printf("Launching GLUT mainloop\n");
    glutMainLoop();
	
	setIntegrator(NULL);

    return EXIT_SUCCESS;
}
