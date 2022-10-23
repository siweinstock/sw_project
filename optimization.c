#include <stdlib.h>
#include <float.h>

#include "spmat.h"
#include "optimization.h"
#include "power_iteration.h"
#include "errors.h"
#define EPSILON 0.00001

/* Optimizes the given division s according to the algorithm */
void optimize(Bmat* B, double* s){
	int i, j, n = B->dim, max_score_ind, max_imp_ind;
	double delta_Q, max_score, max_imp, improve;
	int* indices = (int*)malloc(n * sizeof(int));
	double* score = (double*)malloc(n * sizeof(double));
	double* B_row = (double*)malloc(n * sizeof(double));
	int* vertices_moved = (int*)malloc(n * sizeof(int)); /* This array indicates which vertices were
				    	   	   	   	   	            	  * already moved in this iteration */
	if (indices == NULL || score == NULL || B_row == NULL || vertices_moved == NULL) {
		error_msg(ERR_MEM_ALLOC_FAILED);
	}

	/* This while loop goes until we get to a stage where mod==max_external_state->mod
	 * which means that there was no improvement in this iteration
	 * hence this is the optimal division */
	do{
		for (j=0; j<n; j++){
			vertices_moved[j] = 0; /* Set vertices_moved to be a zeros array */
		}

		/*matrix_vector_mult(B, s, x);*/ /* x = Bs */
		max_score = -DBL_MAX;
		for (j=0; j<n; j++){
			score[j] = -2*(s[j]*compute_Brow_mult_s(B, j, B_row, s)
					+ (double)B->k[j]*(double)B->k[j]/M);
		}

		max_imp = -DBL_MAX;
		improve = 0;
		/* In this loop we do n iterations, one for each vertex,
		 * each iteration we move the vertex that gives us maximum modularity to the other group */
		for (i=0; i<n; i++){
			max_score = -DBL_MAX;

			for (j=0; j<n; j++){

				if (vertices_moved[j]==0){
					if (score[j]>max_score){
						max_score = score[j];
						max_score_ind = j;
					}
				}
			}

			move(s, max_score_ind); /* Moving the vertex to the other group */
			vertices_moved[max_score_ind] = 1; /* Mark the vertex as "moved" */
			compute_B_row(B, max_score_ind, B_row);
			for (j=0; j<n; j++){

				if (j==max_score_ind){
					score[j] = -score[j];
				}else{
					score[j] -= 4 * s[j] * s[max_score_ind] * B_row[j];
				}

			}
			indices[i] = max_score_ind;
			improve += max_score;

			if (improve>max_imp){
				max_imp = improve;
				max_imp_ind = i;
			}
		}

		printf("%f\n", max_imp);
		for (i=n-1; i>max_imp_ind; i--){ /* Set s to be the optimal division */
			move(s, indices[i]);
		}
		if (max_imp_ind == n-1){
			delta_Q = 0;
		}else{
			delta_Q = max_imp;
		}

	}while (delta_Q > EPSILON);
	free(B_row); free(indices);
	free(score); free(vertices_moved);
}


/* Moves the vertex j to the other group */
void move(double* s, int j){
	s[j] *= -1;
}



void compute_B_row(Bmat* B, int i, double* B_row){
	int j;
	B->A->get_row(B->A, i, B_row);
	for (j=0; j<B->A->n; j++){
		B_row[j] -= (double)B->k[i]*(double)B->k[j]/M;
	}
}

double compute_Brow_mult_s(Bmat* B, int i, double* B_row, double* s){
	int j;
	double x = 0;
	B->A->get_row(B->A, i, B_row);
	for (j=0; j<B->A->n; j++){
		B_row[j] -= (double)B->k[i]*(double)B->k[j]/M;
		x += s[j]*B_row[j];
	}
	return x;
}
