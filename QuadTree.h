//
//  QuadTree.h
//  GravSim
//
//  Created by Thomas Dickerson on 4/15/13.
//
//

#ifndef GravSim_QuadTree_h
#define GravSim_QuadTree_h
#include "Constants.h"

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

typedef struct{Vector lowX; Vector highX; } AABB;
typedef struct _qtNode{ AABB bounds; Body * contents; /* contents may be virtual */ struct _qtNode * parent; struct _qtNode * children [2] [2]; } QuadTree;



AABB makeNullBB();
int intersectsAABB(AABB * bb, Vector * v);

Vector dim(AABB * bb);

Vector mid(AABB * bb);

void redimAABB(AABB * bb, AABB * fixed, Vector * mag);
void updateAABB(AABB * bb, Vector * v);

Vector barycenter(Body * b1, Body * b2);

Vector corner(Vector *v1, Vector *v2, int quad);
void updateMasses(QuadTree * node, Body * b);

Body* getDummyBody();

void returnDummyBody(QuadTree * qt);
QuadTree * freeTree(QuadTree * node);

QuadTree * initEmptyQuad();
QuadTree * initInQuad(QuadTree * node, int quad);
void initChildren(QuadTree * node);
void splitNode(QuadTree * node, Body * b);
void insertBody(QuadTree * node, Body * b);

#endif
