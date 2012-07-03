#include <stdio.h>
#include <math.h>
//#include <stdlib.h>

#include <gsl/gsl_matrix.h>

#include "lib_grover_simulation.h"

int main(int argc, char **argv)
{
	/* number of bosons : M
	 * number of sites  : Np
	 * computing from Np_start to Np_stop (arguments from the command line)
	 * The GNU Scientific Library is used for operations on matrices
	 * rows_iterator and columns_iterator and used to run through the matrix H.
	 */

	int M=atoi(argv[3]),
		Np,
		Np_start=atoi(argv[1]),
		Np_stop=atoi(argv[2]),
		rows_iterator,
		columns_iterator,
		combfactor_iterator,
		binomial_iterator;


	int *statens;
	double *combfactor,
		   combfactor_var,
		   combfactor_bin,
		   combfactor_product;

	for(Np=Np_start;Np<=Np_stop;Np++) {
		printf("Np = %d\n",Np);
		// Get maxtrix size
		size_t matrix_size = compute_matrix_size(Np+1,M);
		// Allocate matrix filled of 0s
		gsl_matrix *H = gsl_matrix_calloc(matrix_size,matrix_size);
		// Allocate combfactor array
		combfactor = (double*)malloc(matrix_size*sizeof(double));
		// Iteration over rows of the matrix
		for(rows_iterator=0;rows_iterator<matrix_size;rows_iterator++) {
			statens = compute_state_list(rows_iterator,Np+1,M);
			COMPUTE_COMBFACTOR(Np,statens,M,combfactor[rows_iterator]);
			printf("combfactor[%d] = %.10f\n",rows_iterator,combfactor[rows_iterator]);
			/*
			 *
			 *
			 */
			for(columns_iterator=0;columns_iterator<matrix_size;columns_iterator++) {
				/*
				 *
				 *
				 */

				//free(statens)
				//free(combfactor)
			}
		}

		//gsl_matrix_free(H);
	}
	return 0;
}

