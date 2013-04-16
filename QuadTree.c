//
//  QuadTree.c
//  GravSim
//
//  Created by Thomas Dickerson on 4/15/13.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "QuadTree.h"

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
		if(x > h || isnan(h)) setVField(&(bb->highX), i, x + sqEP);
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
		ob->p = sum(ob->p,b->p);
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
	if(node != NULL)
	{
		returnDummyBody(node);
		FREE_CHILDREN(node);
		free(node);
	}
	return NULL;
}

QuadTree * initEmptyQuad()
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
	QuadTree * ret = initEmptyQuad();
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
	ret->bounds = (AABB){.lowX = lx, .highX = hx};
	return ret;
	
	
}

void initChildren(QuadTree * node)
{
	MAKE_CHILDREN(node);
	
}

void insertBody(QuadTree*,Body*);
void splitNode(QuadTree * node, Body * b)
{
	Body * oldb = node->contents;
	initChildren(node);
	node->contents = getDummyBody();
	Vector nmid = mid(&(node->bounds));
	
	insertBody(node->children[(oldb->x.e1 > nmid.e1)][(oldb->x.e2 > nmid.e2)], oldb);
	insertBody(node->children[(b->x.e1 > nmid.e1)][(b->x.e2 > nmid.e2)], b);
}

void insertBody(QuadTree * node, Body * b)
{
	if(node->contents == NULL){
		node->contents = b;
		if(node->parent != NULL) updateMasses(node->parent, b);
	}
	else if(HAS_NO_CHILDREN(node)) {
		// Nobody here but us chickens
		splitNode(node, b);
	} else{
		Vector nmid = mid(&(node->bounds));
		insertBody(node->children	[(b->x.e1 > nmid.e1)]
				   [(b->x.e2 > nmid.e2)], b);
	}
}