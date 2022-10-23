#include "modularity.h"
#include "network.h"

/* Removes one set from S and returns it */
Community* pop(Network* S){
	Community* s = S->head;
	S->head = s->next;
	S->n--;
	return s;
}

/* Inserts s to S */
void insert(Network* S, Community* s){
	s->next = S->head;
	S->head = s;
	S->n++;
}
