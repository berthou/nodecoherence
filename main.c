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
	 * matelem_iterator is used as iterator to compute the diagonal-specific product 
	 */

	int M=atoi(argv[3]),
		Np,
		Np_start=atoi(argv[1]),
		Np_stop=atoi(argv[2]),
		rows_iterator,
		columns_iterator,
		combfactor_iterator,
		binomial_iterator,
		matelem_iterator,
		matelem_sum_iterator;

	int *statens,
		*statensp;

	size_t matrix_size;

	double *combfactor,
		   combfactor_var,
		   combfactor_bin,
		   combfactor_product,
		   matelem_var,
		   matelem_diagonal_var,
		   matelem_product,
		   matelem_sum,
		   *szmatelem;

	for(Np=Np_start;Np<=Np_stop;Np++) {
		printf("Np = %d\n",Np);
		// Get maxtrix size
		matrix_size = compute_matrix_size(Np+1,M);
		// Allocate matrix
		gsl_matrix *H = gsl_matrix_alloc(matrix_size,matrix_size);
		// Allocate combfactor and szmatelem arrays
		combfactor = (double*)malloc(matrix_size*sizeof(double));
		szmatelem  = (double*)malloc(matrix_size*sizeof(double));
		// Iteration over rows of the matrix
		for(rows_iterator=0;rows_iterator<matrix_size;rows_iterator++) {
			statens = compute_state_list(rows_iterator,Np+1,M);
			COMPUTE_COMBFACTOR(Np,statens,M,combfactor[rows_iterator]);
			szmatelem[rows_iterator]=(2.0*statens[0]-Np)/Np;	

			// Iteration over columns of the matrix
			for(columns_iterator=0;columns_iterator<matrix_size;columns_iterator++) {
				statensp = compute_state_list(columns_iterator,Np+1,M);
				matelem_var=0.0;
				matelem_diagonal_var=1.0;

				// If we are on the diagonal of H :
				if (rows_iterator == columns_iterator) {
					for(matelem_iterator=0;matelem_iterator<M;matelem_iterator++) {
						matelem_diagonal_var*=(double)statens[matelem_iterator]/(double)Np;
					}
					//printf("%f ",matelem_var/*gsl_matrix_get(H,rows_iterator,rows_iterator*/);
				}
				// For every element of the matrix :
				matelem_product=1.0;
				for(matelem_iterator=0;matelem_iterator<M;matelem_iterator++) {
					matelem_sum=0.0;
					for(matelem_sum_iterator=0;matelem_sum_iterator<=Np;matelem_sum_iterator++) {
						matelem_sum+=(((double)matelem_sum_iterator/(double)Np)*inner_product((double)statens[matelem_iterator],(double)matelem_sum_iterator,(double)Np)*inner_product((double)statensp[matelem_iterator],(double)matelem_sum_iterator,(double)Np));
					}
					matelem_product*=matelem_sum;
				}
				// If we are on the diagonal of H (for affectation) :
				if (rows_iterator == columns_iterator)
					gsl_matrix_set(H,rows_iterator,columns_iterator,pow(Np,2)*(matelem_diagonal_var+matelem_product));
				else
					gsl_matrix_set(H,rows_iterator,columns_iterator,pow(Np,2)*(matelem_product));


				//free(statens)
				//free(combfactor)
			}
		}
		print_matrix(H,matrix_size);
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

