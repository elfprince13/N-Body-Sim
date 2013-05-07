//
//  QTNodeList.c
//  GravSim
//
//  Created by Thomas Dickerson on 5/6/13.
//
//

#include <stdlib.h>
#include <stdio.H>
#include "QTNodeList.h"


void prepend(QTList * l, QuadTree * qt){
	QTListNode * qtln = (QTListNode *)malloc(sizeof(QTListNode));
	qtln->prev = NULL;
	qtln->next = NULL;
	qtln->qt = qt;
	if(l->head != NULL){
		l->head->prev = qtln;
		qtln->next = l->head;
		l->head = qtln;
	} else{
		l->head = qtln;
		l->tail = qtln;
	}
}

void append(QTList * l, QuadTree * qt){
	QTListNode * qtln = (QTListNode *)malloc(sizeof(QTListNode));
	qtln->prev = NULL;
	qtln->next = NULL;
	qtln->qt = qt;
	if(l->tail != NULL){
		l->tail->next = qtln;
		qtln->prev = l->tail;
		l->tail = qtln;
	} else{
		l->head = qtln;
		l->tail = qtln;
	}
}

QuadTree* pophead(QTList * l){
	QuadTree* ret = NULL;
	QTListNode *nhead;
	if(l->head != NULL){
		ret = l->head->qt;
		nhead = l->head->next;
		l->head->next = NULL;
		free(l->head);
		l->head = nhead;
		if(nhead == NULL){
			l->tail = NULL;
		}
	}
	return ret;
}

QuadTree* poptail(QTList * l){
	QuadTree* ret = NULL;
	QTListNode *ntail;
	if(l->tail != NULL){
		ret = l->tail->qt;
		ntail = l->head->prev;
		l->tail->prev = NULL;
		free(l->tail);
		l->tail = ntail;
		if(ntail == NULL){
			l->head = NULL;
		}
	}
	return ret;
}

void freeQTList(QTList * l){
	if(l!=NULL){
		while(!isEmptyQTList(l)){
			pophead(l);
		}
		free(l);
	}
}

int isEmptyQTList(QTList *l){	return (l->head == NULL || l->tail == NULL);	}