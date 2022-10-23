#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "power_iteration.h"
#include "errors.h"

#define EPSILON 0.00001

double dot_prod_id(const int* a, const double* b, int n) {
    int i = 0;
    double prod = 0;
    while (i < n) {
        prod += *(a+i) * *(b+i);
        i++;
    }
    return prod;
}

double dot_prod_dd(const double* a, const double* b, int n) {
    int i = 0;
    double prod = 0;
    while (i < n) {
        prod += *(a+i) * *(b+i);
        i++;
    }
    return prod;
}

/*
 * Perform a double vector multiplication by a scalar.
 *
 * */
void scale_double_vector(const double* vec, double scalar, double* result, int n) {
    int i;

    for (i=0; i < n; i++) {
        *(result+i) = *(vec+i) * scalar;
    }
}

/*
 * Perform an int vector multiplication by a scalar.
 *
 * */
void scale_int_vector(const int* vec, double scalar, double* result, int n) {
    int i;

    for (i=0; i < n; i++) {
        *(result+i) = *(vec+i) * scalar;
    }
}


/*
 * Perform an action between the corresponding components of two vectors on n dimensions.
 * Used for vector summation and subtraction and also for the multiplication of a diagonal matrix
 * (which can be represented as a vector of the diagonal) with a vector.
 *
 * Arguments:
 * a, b - vectors with n conponents
 * op - the operation to be performed ('+', '-', '*')
 * n - the vectors' dimension
 *
 * Returns:
 * the result vector for the matching operation
 *
 * */
void vec_componentwise_op(const double* a, const double* b, double* result, char op, int n) {
    int i;

    switch(op) {
        case '+':
            for (i=0; i<n; i++) {
                *(result+i) = *(a+i) + *(b+i);
             }
            break;
        case '-':
            for (i=0; i<n; i++) {
                *(result+i) = *(a+i) - *(b+i);
            }
            break;
        case '*':
            for (i=0; i<n; i++) {
                *(result+i) = *(a+i) * *(b+i);
            }
            break;
        default:
            break;
    }
}


/*
 * Initialize the first eigenvector estimation randomly.
 *
 * Arguments:
 * vec - the eigenvector storage variable
 * n - the vector's dimension
 *
 * Returns:
 * None. (the argument vec is populated with the initial guess)
 *
 * */
void init_vec(double* vec, int n) {
    int i = 0;

    srand(time(NULL));

    while (i < n) {
        *(vec + i) = rand();
        i++;
    }

}


double get_norm_squared(const double* vec, int n) {
    int i;
    double norm_squared = 0;

    for (i=0; i < n; i++) {
        norm_squared += *(vec+i) * *(vec+i);
    }

    return norm_squared;
}


void normalize(double* vec, int n) {
    int i;
    double norm_squared = 0;

    /* Replace with get norm sq */
    for (i=0; i < n; i++) {
        norm_squared += *(vec+i) * *(vec+i);
    }
    /* end replace */

    if (norm_squared == 0) {
        error_msg(ERR_DIVISION_BY_ZERO);
    }

    for (i=0; i < n; i++) {
        *(vec+i) /= sqrt(norm_squared);
    }

}


/*
 * Generates the next guess vector
 *
 * */
void matrix_vector_mult(Bmat* B, double* curr_guess, double* next_vec) {
    double* tmp = (double*) malloc(B->dim * sizeof(double));
    double dot_prod_res;

    if (tmp == NULL) {
        error_msg(ERR_MEM_ALLOC_FAILED);
    }

    B->A->mult(B->A, curr_guess, next_vec);  /* A[g]*x = next_vec */

    dot_prod_res = dot_prod_id(B->k, curr_guess, B->dim);  /* k[g]T*x */

    /* M != 0 guaranteed */
    scale_int_vector(B->k, dot_prod_res / M, tmp, B->dim); /* tmp <- ((k[g]T*x)/M)k[g] */

    vec_componentwise_op(next_vec, tmp, next_vec,'-', B->dim);  /* A[g]*x - (k[g]T*x)k[g] */

    vec_componentwise_op(B->D, curr_guess, tmp, '*', B->dim); /*  tmp <- D*x */

    /*  A[g]*x - ((k[g]T*x)/M)k[g] - D*x */
    vec_componentwise_op(next_vec, tmp, next_vec, '-', B->dim);

    scale_double_vector(curr_guess, B->norm, tmp, B->dim); /* ||B[g]||1*x */

    /*  A[g]*x - ((k[g]T*x)/M)k[g] - D*x + ||B[g]||1*x */
    vec_componentwise_op(next_vec, tmp, next_vec, '+', B->dim);

    free(tmp);

}


/*
 * Check whether the change in eigenvector estimation is small enough.
 *
 * Arguments:
 * vec - previous eigenvector estimation
 * next_vec - current eigenvector estimation
 * n - vectors' dimension
 *
 * Returns:
 * 0 - if ta least one component differ by EPSILON.
 * 1 - if all components differ by less than EPSILON from each other.
 *
 * */
int close_enough(const double* vec, const double* next_vec, int n) {
    int i;

    for (i=0; i < n; i++) {
        if (fabs(*(vec + i) - *(next_vec + i)) >= EPSILON) {
            return 0;
        }
    }

    return 1;
}


/*
 * Calculate the eigenvalue estimation according to the formula:
 *      eigenval = b_(k-1) * b_k / b_(k-1) * b_(k-1)
 *
 * Arguments:
 * prev_vec - previous eigenvector estimation
 * vec - current eigenvector estimation
 * n - vectors' dimension
 *
 * Returns:
 * eigenvalue estimation
 *
 * */
double get_eigenvalue(double* prev_vec, double* vec, int n) {
    double numerator = dot_prod_dd(prev_vec, vec, n);
    double denominator = dot_prod_dd(prev_vec, prev_vec, n);

    if (denominator == 0) {
        error_msg(ERR_DIVISION_BY_ZERO);
    }

    return numerator / denominator;
}

/*
 * Compute the leading eigenpair for a modularity matrix.
 *
 * Arguments:
 * B - modularity matrix
 *
 * Returns:
 * the leading eigenvalue and eigenvector packed a an EigenPair struct
 *
 * */
EigenPair* find_leading_eigenpair(Bmat* B) {
    double *vec;
    double* eigenvec;
    double* tmp; /* temporary vector placeholder */
    int iters = 0; /* count number of iterations */
    double eigenvalue;
    int n = B->dim;
    EigenPair* result = (EigenPair*) malloc(sizeof(EigenPair));
    eigenvec = (double*) malloc(n * sizeof(double));
    vec = (double*) malloc(n*sizeof(double));

    if (result == NULL || eigenvec == NULL || vec == NULL) {
        error_msg(ERR_MEM_ALLOC_FAILED);
    }

    init_vec(eigenvec, n); /* Guess initial vector */
    matrix_vector_mult(B, eigenvec, vec);
    eigenvalue = get_eigenvalue(vec, eigenvec, n);
    normalize(vec, n);

    /* Keep guessing new vectors until change in each component is < EPSILON */
    while (!close_enough(eigenvec, vec, n)) {
        tmp = eigenvec;
        eigenvec = vec;
        vec = tmp;

        matrix_vector_mult(B, eigenvec, vec);
        eigenvalue = get_eigenvalue(eigenvec, vec, n) - B->norm;
        normalize(vec, n);

        iters++;
    }

    /* printf("%d iterations\n", iters); */

    result->eigenvalue = eigenvalue;
    result->eigenvector = vec;

    free(eigenvec);

    return result;
}


int free_ep(EigenPair* ep) {
    free(ep->eigenvector);
    free(ep);
    return 0;
}
