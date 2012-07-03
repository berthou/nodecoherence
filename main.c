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
	 * rows_iterator and columns_iterator and used to run through the matrix H
	 * binomial_iterator and combfactor_iterator are respectively used in the BINOMIAL_COEFF and COMPUTE_COMBFACTOR macros as iterators
	 * combfactor_var, combfactor_bin, combfactor_product are variables of the COMPUTE_COMBFACTOR macro
	 * diagonal_iterator is used as iterator to compute the diagonal-specific product 
	 */

	int M=atoi(argv[3]),
		Np,
		Np_start=atoi(argv[1]),
		Np_stop=atoi(argv[2]),
		rows_iterator,
		columns_iterator,
		combfactor_iterator,
		binomial_iterator,
		diagonal_iterator;


	int *statens,
		*statensp;
	
	double *combfactor,
		   combfactor_var,
		   combfactor_bin,
		   combfactor_product,
		   diagonal_product,
		   *szmatelem;

	for(Np=Np_start;Np<=Np_stop;Np++) {
		printf("Np = %d\n",Np);
		// Get maxtrix size
		size_t matrix_size = compute_matrix_size(Np+1,M);
		// Allocate matrix filled of 0s
		gsl_matrix *H = gsl_matrix_calloc(matrix_size,matrix_size);
		// Allocate combfactor array
		combfactor = (double*)malloc(matrix_size*sizeof(double));
		szmatelem  = (double*)malloc(matrix_size*sizeof(double));
		// Iteration over rows of the matrix
		for(rows_iterator=0;rows_iterator<matrix_size;rows_iterator++) {
			statens = compute_state_list(rows_iterator,Np+1,M);
			COMPUTE_COMBFACTOR(Np,statens,M,combfactor[rows_iterator]);
			szmatelem[rows_iterator]=(2.0*statens[0]-Np)/Np;	

			for(columns_iterator=0;columns_iterator<matrix_size;columns_iterator++) {
				statensp = compute_state_list(columns_iterator,Np+1,M);
				
				// If we are on the diagonal of H :
				if (rows_iterator == columns_iterator) {
						diagonal_product=1.0;
					for(diagonal_iterator=0;diagonal_iterator<M;diagonal_iterator++) {
						diagonal_product*=(float)statens[diagonal_iterator]/(float)Np;
					}
					gsl_matrix_set(H,rows_iterator,columns_iterator,diagonal_product);
					//printf("%f ",diagonal_product/*gsl_matrix_get(H,rows_iterator,rows_iterator*/);
				}
				// For every element of the matrix :
					//gsl_matrix_set(H,rows_iterator,columns_iterator,

				//free(statens)
				//free(combfactor)
			}
		}
		/*
		int i;
		for(i=0;i<M;i++) {
			printf("statensp[%d] = %d\n",i,statensp[i]);
			printf("statens[%d] = %d\n",i,statens[i]);
		}
		*/
		//gsl_matrix_free(H);
	}
	return 0;
}

