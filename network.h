#ifndef NETWORK_H_
#define NETWORK_H_

#include "modularity.h"

/* A struct for the sets P and O */
typedef struct {
	Community* head;
	int n;
}Network;

Community* pop(Network* S);

void insert(Network* S, Community* s);

#endif /* NETWORK_H_ */
