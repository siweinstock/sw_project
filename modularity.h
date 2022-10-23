/*
 * Modularity.h
 *
 * This module holds the functions necessary in order to build modularity matrices B
 * as well as calculating their modularity Q.
 *
 *
 * */

#ifndef ON_MODULARITY_H
#define ON_MODULARITY_H

#include <stdio.h>
#include "spmat.h"

/* Global variables for algorithm's functionality */
int M; /* Sum of entire graph's node degrees */
int* k; /* The node degree vector for the entire graph */
int dim; /* The total number of nodes in the graph (= matrices and vectors degree) */

typedef struct state {
    int row;
    int neighbor_count;
    int* neighbor_list;
    double* row_data;
} State;

State* init_state(int community_size);
void reset_state(State* state);
void free_state(State* state);

/*
 * Community - a structure to represent a group of nodes within a graph
 *
 * Fields:
 * size - number of nodes in community
 * nodes - list of nodes in the community
 * next - a pointer to the next community (in queue for processing)
 *
 * */
typedef struct _Community {
    int size;
    int* nodes;
    struct _Community* next;
} Community;


void free_community(Community* g);


/*
 * Bmat - a structure to represent a modularity matrix
 *
 * Fields:
 * dim - matrix dimension (= number of nodes in community represented by this matrix)
 * A - sparse adjacency matrix between nodes in the community
 * D - diagonal matrix that holds the sum of current row in the B matrix
 * k - community node degree vector (built from global k vector)
 *
 * */
typedef struct _B {
    int dim;
    spmat* A;
    double* D;
    int* k;
    double norm;
} Bmat;


void free_B(Bmat* B);


/*
 * Initializes the global variables in this module based on data from input file
 *
 * Arguments:
 * f_in - input graph binary file
 *
 * Returns:
 * None. (Global variables' values are set)
 *
 * */
void initialize(FILE* f_in);

Community* trivial_community();

/*
 * Build a modularity matrix for the current community.
 *
 * Arguements:
 * f_in - the input graph binary file.
 * k - the entire graph's k-vector
 * g - the current community for which we build the community matrix.
 *
 * Returns:
 * B - the modularity matrix for the community g.
 *
 * */
Bmat* community_modularity_matrix(FILE* f_in, Community* g, State* state);

/*
 * Generate the vector s which resresent a division of the current community into 2 according to the modularity matrix
 * leading eigenvector.
 *
 * Arguments:
 * eigenvector - the modularity matrix B's leading eigenvector.
 * n - the number of components in the eigenvector.
 *
 * Returns:
 * a vector of {1,-1} which specifies to which group in the division the corresponding node belongs to.
 *
 * */
void compute_s(const double* eigenvector, int n, double* s);


Community* debug_community();
Community* debug_community2();
Community* debug_community3();

/*
 * Compute the modularity of the community according to the formula:
 *      Q = s^T * B * s
 *
 * The calculation is done separately for each component of B in order to improve time complexity.
 *
 * Arguments:
 * B - the community's modularity matrix.
 * s - the community division into 2 vector.
 * g - the community.
 *
 * Returns:
 * Q - the community's modularity
 *
 * */
double compute_modularity(Bmat* B, const double* s, Community* g);

/*
 *
 *
 * */
void divide(Bmat* B, Community* g, Community** division);

void free_division(Community** division);

#endif /*ON_MODULARITY_H*/
