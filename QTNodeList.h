//
//  QTNodeList.h
//  GravSim
//
//  Created by Thomas Dickerson on 5/6/13.
//
//

#ifndef GravSim_QTNodeList_h
#define GravSim_QTNodeList_h

#include "QuadTree.h"

typedef struct _qtLN{struct _qtLN * prev;struct _qtLN * next; QuadTree * qt; } QTListNode;
typedef struct _qtL{QTListNode * head; QTListNode * tail;} QTList;

void prepend(QTList * l, QuadTree * qt);
void append(QTList * l, QuadTree * qt);

QuadTree* pophead(QTList * l);
QuadTree* poptail(QTList * l);
void freeQTList(QTList * l);
int isEmptyQTList(QTList *l);


#endif
