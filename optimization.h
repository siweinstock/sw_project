/*
 * optimization.h
 *
 *  Created on: 30 ����� 2020
 *      Author: alonb
 */

#ifndef OPTIMIZATION_H_
#define OPTIMIZATION_H_

#include "modularity.h"


void move(double* s, int j);/* Moves the vertex j to the other group */

/* Optimizes the given division s according to the algorithm */
void optimize(Bmat* B, double* s);

/* Computes the row j of B, Bj */
void compute_B_row(Bmat* B, int j, double* Bj);

double compute_Brow_mult_s(Bmat* B, int i, double* B_row, double* s);

#endif /* OPTIMIZATION_H_ */
