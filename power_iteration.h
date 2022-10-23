#ifndef ON_POWER_ITERATION_H
#define ON_POWER_ITERATION_H
#include "modularity.h"

typedef struct _eigenpair {
    double eigenvalue;
    double* eigenvector;
} EigenPair;

int free_ep(EigenPair* ep);

void init_vec(double* vec, int n);
EigenPair* find_leading_eigenpair(Bmat* B);


double dot_prod_id(const int* a, const double* b, int n);
void scale_int_vector(const int* vec, double scalar, double* result, int n);
void vec_componentwise_op(const double* a, const double* b, double* result, char op, int n);
void scale_double_vector(const double* vec, double scalar, double* result, int dim);
void normalize(double* vec, int n);

void matrix_vector_mult(Bmat* B, double* curr_guess, double* next_vec);
double dot_prod_dd(const double* a, const double* b, int n);
double get_norm_squared(const double* vec, int n);

#endif /*ON_POWER_ITERATION_H*/
