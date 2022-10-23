#include <stdio.h>
#include <stdlib.h>

#include "spmat.h"
#include "errors.h"

typedef struct Cell {
    int row;
    int col;
    double val;
    struct Cell* next;
} Cell;


Cell* new_cell(int row, int col, double val) {
    Cell* cell = (Cell*) malloc(sizeof(Cell));
    cell->next = NULL;
    cell->row = row;
    cell->col = col;
    cell->val = val;

    if (cell == NULL) {
        error_msg(ERR_MEM_ALLOC_FAILED);
    }

    return cell;
}

void print_list(Cell* head) {
    while (head != NULL) {
        printf("[%d, %d]%.2f, ", head->row, head->col, head->val);
        head = head->next;
    }
    printf("\n");
}


void cleanup(Cell* head) {
    if (head != NULL) {
        cleanup(head->next);
        free(head);
    }
}

/*****************************************************************************************************************/

void sp_print(struct _spmat *A) {
    int dim = A->n;
    int i;

    for (i=0; i<dim; i++) {
        print_list(*((Cell**)A->private + i));
    }

}

void ll_add_row(struct _spmat *A, const double *row, int i) {
    int dim = A->n;
    int j; /* col num */
    int is_first = 1;
    Cell *first, *last, *tmp;

    for (j=0; j<dim; j++) {
        if(*(row+j) != 0) {
            if (is_first) {
                first = new_cell(i, j, *(row+j));
                last = first;
                is_first = 0;
            }
            else {
                tmp = new_cell(i, j, *(row+j));
                last->next = tmp;
                last = tmp;
            }
        }
    }

    /* If current row contains only zeros */
    if (is_first) {
        *(((Cell**)A->private) + i) = NULL;
    }
    /* If current row contains only zeros */
    else {
        *(((Cell**)A->private) + i) = first;
    }
}

void ll_mult(const struct _spmat *A, const double *v, double *result) {
    int dim = A->n;
    int i, row, col;
    double val;
    int count;

    Cell** lists = (Cell**)A->private;
    Cell* cell;

    for(i=0; i<dim; i++) { /* iterate matrix rows */
        cell = (*(lists+i));
        *(result+i) = 0;
        count = 0;

        while(cell != NULL) { /* multiply non-zero values */
            row = cell->row;
            col = cell->col;
            val = cell->val;

            *(result+row) += val * *(v+col);

            cell = cell->next;
            count++;
        }
    }



}

void ll_free(struct _spmat *A) {
    int dim = A->n;
    int i;
    Cell** lists = (Cell**) A->private;

    /* for each list in array free all cells */
    for (i=0; i<dim; i++) {
        cleanup(*(lists+i));
    }

    free(lists); /* free array */
    free(A); /* free spmat */

}

void ll_get_row(const struct _spmat *A, int i, double* row){
	int j;
	Cell** lists = (Cell**) A->private;
	Cell* current_cell = (Cell*) lists[i];



	for (j=0; j<A->n; j++){
		row[j]=0;
	}


	while (current_cell != NULL){
		row[current_cell->col] = current_cell->val;
		current_cell = current_cell->next;
	}

}

spmat* spmat_allocate_list(int n) {
    spmat* p_spmat = (spmat*) malloc(sizeof(spmat));
    Cell* p_list = malloc(n*sizeof(Cell)); /* initialize an array of n Cell pointers */

    if (p_spmat == NULL || p_list == NULL) {
        error_msg(ERR_MEM_ALLOC_FAILED);
    }

    p_spmat->n = n;
    p_spmat->private = p_list;
    p_spmat->add_row = ll_add_row;
    p_spmat->mult = ll_mult;
    p_spmat->get_row = ll_get_row;
    p_spmat->free = ll_free;
    p_spmat->print = sp_print;

    return p_spmat;
}/**/

/*****************************************************************************************************************/

/* Defining an inner structure for matrices with arrays implementation */
typedef struct{

    double* values; /* An array for the non-zero values */

    int* colind; /* An array for the non-zero values' column index */

    int* rowptr; /* An array for first non-zero value index after each row */

    int counter; /* For add_row */

}array_mat;





/* Adds row i to the matrix. Called before any other call,
 * exactly n times in order (i = 0 to n-1), Implemented for array spmat */
void array_add_row(struct _spmat *A, const double *row, int i){
    int j, n;


    n = A->n; /* Extracting n from A */


    *((*(array_mat*) (A->private)).rowptr+i) = (*(array_mat*) (A->private)).counter;

    /* For each element in row check if non-zero and if so,
     * put the non-zero value in values array and the index in colind array. */
    for (j=0; j<n; j++){
        if (*(row+j) != 0.0){
            *((*(array_mat*) (A->private)).values+(*(array_mat*) (A->private)).counter) = *(row+j);
            *((*(array_mat*) (A->private)).colind+(*(array_mat*) (A->private)).counter) = j;
            (*(array_mat*) (A->private)).counter++;
        }
    }


}


/* Frees all resources used by A , Implemented for array spmat */
void array_free(struct _spmat *A){
    array_mat mat = *(array_mat*) (A->private);
    free(mat.values); free(mat.colind); free(mat.rowptr);
    free(A->private);
    free(A);
}


/* Multiplies matrix A by vector v, into result (result is pre-allocated),
 *  Implemented for array spmat */
void array_mult(const struct _spmat *A, const double *v, double *result){
    int cnt, i, j, n, m;
    double sum;


    n = A->n; /* Extracting n from A */


    cnt = 0;
    /* Calculates each element and puts it in result */
    for (i=0; i<n; i++, cnt += m){
        sum = 0;
        m = (*(array_mat*) (A->private)).rowptr[i+1] - (*(array_mat*) (A->private)).rowptr[i];
        for (j=cnt; j<cnt+m; j++){
            sum += v[(*(array_mat*) (A->private)).colind[j]] * ((*(array_mat*) (A->private)).values[j]);
        }

        result[i] = sum;
    }
}


/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat* spmat_allocate_array(int n, int nnz) {

    array_mat* mat;
    spmat* spmatptr = (spmat*) malloc(sizeof(spmat));


    mat = (array_mat*) malloc(sizeof(array_mat));

    mat->values = (double*)malloc(nnz * sizeof(double)); /* Allocates a new array for values */

    mat->colind = (int*)malloc(nnz * sizeof(int)); /* Allocates a new array for colind */

    mat->rowptr = (int*)malloc((n + 1) * sizeof(int)); /* Allocates a new array for rowptr */

    mat->counter = 0;
    *(mat->rowptr+n) = nnz; /* The last entry of rowptr is number of non-zero elements */
    spmatptr->n = n; /* n^2 is the size of the matrix */
    spmatptr->private = mat; /* Put private as a pointer to mat */
    spmatptr->add_row = array_add_row;
    spmatptr->mult = array_mult;
    spmatptr->free = array_free;

    return spmatptr;

}

