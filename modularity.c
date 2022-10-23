#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "modularity.h"
#include "power_iteration.h"
#include "optimization.h"
#include "errors.h"


/*typedef struct state {
    int row;
    int neighbor_count;
    int* neighbor_list;
    double* row_data;
} State;*/


State* init_state(int community_size) {
    State* state = (State*) malloc(sizeof(State));

    state->row = 0;
    state->row_data = (double*) malloc(community_size * sizeof(double));
    state->neighbor_list = (int*) malloc(dim * sizeof(int));

    return state;
}

void reset_state(State* state) {
    state->row = 0;
}


void free_state(State* state) {
    free(state->neighbor_list);
    free(state->row_data);
    free(state);
}


/*
 * Build an array of node degree k_i for each node in graph
 *
 * The function reads the 'graph.in' file once and reads only node degree values.
 * The function returns an array of ints of the node degrees.
 *
 * Arguments:
 * f_in - 'graph.in' file
 * dim - matrix dimension
 *
 * Return:
 * k - nodes' degrees vector
 *
 * */
int* build_k_vector(FILE* f_in) {
    int i;
    int* res;

    if ((res = (int*) calloc(dim, sizeof(int))) == NULL) {
        error_msg(ERR_MEM_ALLOC_FAILED);
    }
    if (fseek(f_in, sizeof(int), SEEK_SET)) {
        error_msg(ERR_FSEEK_FAILED);
    }

    for (i=0; i<dim; i++) {
        if (fread(res + i, sizeof(int), 1, f_in) != 1) {
            error_msg(ERR_READ_WRITE_FAILED);
        }
        if (fseek(f_in, *(res + i) * sizeof(int), SEEK_CUR)) {
            error_msg(ERR_FSEEK_FAILED);
        }
    }

    return res;
}

int vec_sum_components(const int* vec, int n) {
    int i, sum = 0;
    for (i=0; i<n; i++) {
        sum+= *(vec + i);
    }
    return sum;
}

/*
 * Initializes global variables
 *
 * */
void initialize(FILE* f_in) {
    if (fread(&dim, sizeof(int), 1, f_in) != 1) {
        error_msg(ERR_READ_WRITE_FAILED);
    }
    k = build_k_vector(f_in);
    M = vec_sum_components(k, dim);

    if (M == 0) { /* no error right now but there will be later on */
        error_msg(ERR_DIVISION_BY_ZERO);
    }
}


void init_B(Bmat* B, Community* g) {
    B->A = spmat_allocate_list(g->size);
    B->D = (double*) malloc(g->size * sizeof(double));
    B->k = (int*) malloc(g->size * sizeof(int));
    B->norm = 0;
    B->dim = g->size;
}


Community* trivial_community() {
    int i;
    int* nodes = (int*) malloc(dim * sizeof(int));
    Community *g = (Community*) malloc(sizeof(Community));

    if (nodes == NULL || g == NULL) {
        error_msg(ERR_MEM_ALLOC_FAILED);
    }

    for (i=0; i<dim; i++) {
        *(nodes+i) = i;
    }

    g->nodes = nodes;
    g->size = dim;
    g->next = NULL;

    return g;
}


/* Runs on rows in group
 * Helps constructing the current row in the modularity's matrix.
 *
 *
 * */
void modularity_row(FILE* f_in, Community* g, Bmat* B, int row, State* state) {
    int group_offset;
    int neighbor_offset = 0;
    int Agi_sum = 0;
    int kgj_sum = 0;
    double curr;
    double row_sum = 0;
    double row_norm1 = 0;

    int* neighbor_list = state->neighbor_list;
    double* row_data = state->row_data;
    int neighbor_count = state->neighbor_count;
    int Ag_row = state->row;

    if (neighbor_list == NULL || row_data == NULL) {
        error_msg(ERR_MEM_ALLOC_FAILED);
    }

    if ((int)fread(neighbor_list, sizeof(int), neighbor_count, f_in) != neighbor_count) {
        error_msg(ERR_READ_WRITE_FAILED);
    }

    for (group_offset=0; group_offset<g->size; group_offset++) { /* O(n) */
        /* can assume less than O(n) because adjacency matrix is sparse */
        while (neighbor_offset < neighbor_count-1 &&
                *(neighbor_list+neighbor_offset) < *(g->nodes+group_offset)) { /* O(1) */
            neighbor_offset++;
        }
        /* if current neighbor is in community */
        if (neighbor_offset < neighbor_count && *(neighbor_list+neighbor_offset) == *(g->nodes+group_offset)) { /* O(1) */
            *(row_data+group_offset) = 1;
            Agi_sum++;
            curr = 1;
        }
        else {
            *(row_data+group_offset) = 0;
            curr = 0;
        }

        curr -= ((double)*(k+row) * *(k+*(g->nodes+group_offset))) / M;
        row_sum += curr;
        kgj_sum += *(k+*(g->nodes+group_offset));

        if (row != *(g->nodes+group_offset)) { /* O(1) */
            row_norm1 += fabs(curr);
        }
    }

    /* M != 0 at this point guaranteed */
    *(B->D+Ag_row) = Agi_sum - (double)(*(k+row) * kgj_sum) / M;
    row_norm1 += (((double) *(k+row) * *(k+row)) / M + row_sum);
    B->A->add_row(B->A, row_data, Ag_row);
    *(B->k+Ag_row) = *(k+row);

    if (row_norm1 > B->norm) {
        B->norm = row_norm1;
    }

}


/*
 * Build the community's modularity matrix.
 * The matrix consists 4 parts:
 * A - the adjacency matrix of the community nodes
 * k - the community node degree vector
 * D - a diagonal matrix where the i'th row is the sum off the modularity matrix B i'th row
 * norm1 - the matrix's norm 1 (sum of column's absolute values)
 *
 * Arguments:
 * f_in - the input graph's binary file.
 * k - the global k vector (for the entire graph).
 * g - the current community.
 *
 * Returns:
 * B - the community's modularity matrix.
 *
 * */
Bmat* community_modularity_matrix(FILE* f_in, Community* g, State* state) {
    /* Modularity matrix components */
    Bmat* B = (Bmat*) malloc(sizeof(Bmat));
/*
    State* state;
*/

    /* Row pointers and data variables */
    int row;
    int node_offset = 0;

    /* create a B object and return it */
    init_B(B, g);
    if (B->A == NULL || B == NULL || B->D == NULL || B->k == NULL) {
        error_msg(ERR_MEM_ALLOC_FAILED);
    }

/*
    state = init_state(g->size);
*/

    if(fseek(f_in, sizeof(int), SEEK_SET)) {
        error_msg(ERR_FSEEK_FAILED);
    }

    for (row=0; row < dim; row++) { /* go over entire graph's rows */
        if ((fread(&state->neighbor_count, sizeof(int), 1, f_in)) != 1) {
            error_msg(ERR_READ_WRITE_FAILED);
        }

        if (node_offset < g->size && row == *(g->nodes + node_offset)) {
            modularity_row(f_in, g, B, row, state);
            state->row++;
            node_offset++;
        }
        else {
            if (fseek(f_in, state->neighbor_count * sizeof(int), SEEK_CUR)) {
                error_msg(ERR_FSEEK_FAILED);
            }
        }
    }

/*
    free_state(state);
*/
    return B;
}


/*
 * Compute the division vector according to the leading eigenvector
 *
 * Arguments:
 * eigenvector - the modularity's matrix leading eigenvector
 * n - the dimension of the eigenvector (= number of nodes in community)
 *
 * Returns:
 * s - the division vector
 *
 * */
void compute_s(const double* eigenvector, int n, double* s) {
    int i;

    if (s == NULL) {
        error_msg(ERR_MEM_ALLOC_FAILED);
    }

    for (i=0; i<n; i++) {
        if (*(eigenvector+i) > 0) {
            *(s+i) = 1;
        }
        else {
            *(s+i) = -1;
        }
    }

}


/*
 * Compute a community division's modularity.
 *
 * Arguments:
 * B - the communities modularity matrix
 * s - the community division vector
 * g - the current community
 *
 * Returns:
 * Q - the community's modularity
 *
 * */
double compute_modularity(Bmat* B, const double* s, Community* g) {
    int i;
    int n = g->size;
    double* sg = (double*) malloc(n * sizeof(double));
    double* kg = (double*) malloc(n * sizeof(double));
    double* sgAg = (double*) malloc(n * sizeof(double));
    double kgsg;
    double Q;

    double* res = (double*) malloc(n*sizeof(double));

    if (sg == NULL || kg == NULL || sgAg == NULL || res == NULL) {
        error_msg(ERR_MEM_ALLOC_FAILED);
    }

    for (i=0; i<n; i++) {
        *(sg+i) = *(s+i);
        *(kg+i) = *(k + *(g->nodes+i));
    }

    for (i=0; i<n; i++) {
        B->A->mult(B->A, sg, sgAg); /* A is symmetric */
    }

    /* M != 0 guaranteed */
    Q = dot_prod_dd(sgAg, sg, n); /* sAs */
    kgsg = dot_prod_dd(sg, kg, n);
    Q -= (kgsg * kgsg / M); /* sAs - sT*(k*kT/M)*s */
    vec_componentwise_op(sg, B->D, res, '*', n);
    Q -= dot_prod_dd(res, sg, n); /* sAs - sT*(k*kT/M)*s - sDs*/

    free(sg);
    free(kg);
    free(sgAg);
    free(res);

    return Q/2;
}


/*
 * Compute a division of a community into 2 given its modularity matrix.
 *
 * Arguments:
 * B - the community's modularity matrix.
 * g - the current community
 * division - an array where the 2 sub-communities go
 *
 * Return:
 * None. (the division is stored at the division argument).
 *
 * */
void divide(Bmat* B, Community* g, Community** division) {
    EigenPair* ep;
    double* sg = (double*) malloc(g->size * sizeof(double));
    double Q;
    int i;
    Community *g1, *g2;
    int  g1_count, g2_count;
    int *g1_nodes, *g2_nodes;
    int g1_p, g2_p;

    ep = find_leading_eigenpair(B);
    compute_s(ep->eigenvector, B->dim, sg);

    free_ep(ep);

    optimize(B, sg);

    Q = compute_modularity(B, sg, g);


    if (Q < 0) {
        division[0] = NULL;
        division[1] = NULL;
    }
    else {
    	g1 = (Community*) malloc(sizeof(Community));
    	g2 = (Community*) malloc(sizeof(Community));

    	if (g1 == NULL || g2 == NULL) {
    	    error_msg(ERR_MEM_ALLOC_FAILED);
    	}

/*
    	optimize(B, sg, g);
*/
        g1_count = g2_count = 0;
        g1_p = g2_p = 0;
        for (i=0; i<B->dim; i++) {
            *(sg + i) == 1 ? g1_count++ : g2_count++ ;
        }

        g1_nodes = (int*) malloc(g1_count * sizeof(int));
        g2_nodes = (int*) malloc(g2_count * sizeof(int));

        if (g1_nodes == NULL || g2_nodes == NULL) {
            error_msg(ERR_MEM_ALLOC_FAILED);
        }


        for (i=0; i<B->dim; i++) {
            if (*(sg + i) == 1) {
                *(g1_nodes+g1_p) = g->nodes[i];
                g1_p++;
            }
            else {
                *(g2_nodes+g2_p) = g->nodes[i];
                g2_p++;
            }
        }

        g1->size = g1_count;
        g1->nodes = g1_nodes;
        g2->size = g2_count;
        g2->nodes = g2_nodes;

        division[0] = g1;
        division[1] = g2;
    }
    free(sg);
}


void free_B(Bmat* B) {
    B->A->free(B->A);
    free(B->D);
    free(B->k);
    free(B);
}


