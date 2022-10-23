#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "modularity.h"
#include "network.h"
#include "errors.h"


int main(int argc, char* argv[]) {

	FILE *f_g;
	FILE *output;
	int i;
	Bmat* B;
	Community *group, *temp;
	Network* P = (Network*)malloc(sizeof(Network));
	Network* O = (Network*)malloc(sizeof(Network));
	Community* trivial_division; /* This represents the group of all
	                              * vertices that is initially placed in P */

	State* state;
	clock_t start, end;

	start = clock();

	if (argc != 3) {
        error_msg(ERR_INVALID_ARG_NUM);
	}

	if ((f_g = fopen(argv[1], "rb")) == NULL) {
	    error_msg(ERR_INVALID_IN_FILE);
	}

    if ((output = fopen(argv[2], "wb")) == NULL) {
        error_msg(ERR_INVALID_IN_FILE);
    }

	initialize(f_g);

	P->head = NULL;
	O->head = NULL;
	P->n = 0;
	O->n = 0;

	trivial_division = trivial_community();

	insert(P, trivial_division);

	state = init_state(dim);

	/* The while loop goes until P is empty, each iteration removes one group g out of P,
	 * then uses algo2 to find the optimal division of g. */
	while (P->head != NULL) {
		Community* g = pop(P); /* Removes g from P */
		Community** division = (Community**) malloc(2 * sizeof(Community*));
		int n1, n2;

		if (division == NULL) {
            error_msg(ERR_MEM_ALLOC_FAILED);
		}

		reset_state(state);
		B = community_modularity_matrix(f_g, g, state);
		divide(B, g, division);



		if (division[0] != NULL) {
            n1 = division[0]->size;
            n2 = division[1]->size;
		}

        /* If one of the groups is empty it means that the optimal division is the trivial division into
         * one group, so we put g in O. Otherwise, we create the groups g1 and g2, put in them the appropriate
         * vertices and check if they are of size 1 */
		if (division[0] == NULL || n1 == 0 || n2 == 0) {
			insert(O, g);
			if (division[0]!=NULL){
			    free(division[0]->nodes);
			    free(division[1]->nodes);
			}
			free(division[0]); free(division[1]);
		}
		else {
			free(g->nodes);
			free(g);

            /* If a group is of size 1 it means it can no longer be divided so we put it in O */
			if (n1 == 1) {
				insert(O, division[0]);
			}
			else {
				insert(P, division[0]);
			}

            if (n2 == 1){
				insert(O, division[1]);
			}
            else {
				insert(P, division[1]);
			}


        }

		free_B(B);
        free(division);

    }

    fclose(f_g);
    free(P);
	free(k);
	free_state(state);

    /* Here we write the groups we found to the output file */
	if (fwrite(&O->n, sizeof(int), 1, output) != 1) {
	    error_msg(ERR_READ_WRITE_FAILED);
	}

    printf("Number of groups = %d\n", O->n);

	group = O->head;

    while (group != NULL) {
		if (fwrite(&(group->size), sizeof(int), 1, output) != 1) {
            error_msg(ERR_READ_WRITE_FAILED);
		}

		printf("Group size = %d\n Group vertices: ", group->size);

		if ((int)fwrite(group->nodes, sizeof(int), group->size, output) != group->size) {
            error_msg(ERR_READ_WRITE_FAILED);
		}

		for (i=0; i<group->size; i++) {
			printf("%d ", group->nodes[i]);
		}
		printf("\n");

		temp = group;
		group = group->next;
		free(temp->nodes);
		free(temp);
	}
	free(O);

	fclose(output);

	end = clock();
	printf("elapsed time: %f\n", (double)(end-start)/CLOCKS_PER_SEC);

	return 0;

}
