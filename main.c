#include <stdio.h>
#include <math.h>

#include "lib_grover_simulation.h"

/* 
 *
 *		!!! check return values
 *
 */

int main(int argc, char **argv)
{
	/* number of bosons : M
	 * number of sites  : Np
	 * computing from Np_start to Np_stop (arguments from the command line)
	 * The GNU Scientific Library is used for operations on matrices
	 * binomial_iterator and combfactor_iterator are respectively used in the BINOMIAL_COEFF and COMPUTE_COMBFACTOR macros as iterators
	 * combfactor_bin, combfactor_product are variables of the COMPUTE_COMBFACTOR macro
	 */

	int M=atoi(argv[3]),
		Np,
		Np_start=atoi(argv[1]),
		Np_stop=atoi(argv[2]),
		matrix_size,
		i,j,n,knx,
		notpoints=40,
		binomial_iterator,
		combfactor_iterator;

	double *statens    = NULL,
		   *statensp   = NULL,
		   *points     = NULL,
		   *combfactor = NULL,
		   *szmatelem  = NULL,
		   *tempvector = NULL,
		   diag_var,
		   matelem,
		   sum,
		   tmax=30,
		   combfactor_bin,
		   combfactor_product;


	gsl_matrix *H,
			   *vectors;

	gsl_eigen_symmv_workspace *eigen_workspace;
	gsl_vector *energies;

#ifdef BENCHMARK
	printf("Starting at : ");
	fflush(stdout);
	system("date");
#endif

	/* Allocate statens and statensp arrays */
	statens  = (double*)malloc(M*sizeof(double));
	statensp = (double*)malloc(M*sizeof(double));

	for(Np=Np_start;Np<=Np_stop;Np++) {
#ifndef BENCHMARK
		printf("Np = %d\n",Np);
#endif
		/* Get matrix size */
		matrix_size=pow(Np+1,M);
		/* Allocate matrix */
		H = gsl_matrix_alloc(matrix_size,matrix_size);
		/* Allocate combfactor, szmatelemarrays */
		combfactor = (double*)realloc(combfactor,matrix_size*sizeof(double));
		szmatelem  = (double*)realloc(szmatelem,matrix_size*sizeof(double));
		/* Iteration over rows of the matrix */
		for(i=0;i<matrix_size;i++)
		{
			compute_state_list(i,Np+1,M,statens);
			COMPUTE_COMBFACTOR(Np,statens,M,combfactor[i]);
			szmatelem[i]=(2.0*statens[0]-(double)Np)/(double)Np;

			/* Iteration over columns of the matrix */
			for(j=0;j<matrix_size;j++)
			{
				compute_state_list(j,Np+1,M,statensp);

				/* If we are on the diagonal of H : */
				if(i==j)
				{
					diag_var=1.0;
					for(n=0;n<M;n++)
					{
						diag_var*=(statens[n]/(double)Np);
					}
				}
				/* For the upper part of the matrix : */
				if(j>=i) 
				{
					matelem=1.0;
					for(n=0;n<M;n++)
					{
						sum=0.0;
						for(knx=0;knx<=Np;knx++)
						{
							sum+=((double)knx/(double)Np)*inner_product(statens[n],(double)knx,(double)Np)*inner_product(statensp[n],(double)knx,(double)Np);
						}
						matelem*=sum;
					}
					if (i==j)
					{
						gsl_matrix_set(H,i,i,pow(Np,2)*(matelem+diag_var));
					}
					else
					{
						gsl_matrix_set(H,i,j,pow(Np,2)*matelem);
						gsl_matrix_set(H,j,i,pow(Np,2)*matelem);
					}
				}
			}
		}
#ifdef DEBUG
		print_matrix(H,matrix_size);
#endif

		/* Allocate workspace needed by the gsl library */
		eigen_workspace = gsl_eigen_symmv_alloc(matrix_size);
		/* Allocate matrix and vector : */
		energies = gsl_vector_alloc(matrix_size);
		vectors = gsl_matrix_alloc(matrix_size,matrix_size);
		/* Get eigenvalues (energies) and eigenvectors (vectors) of H : */
		gsl_eigen_symmv(H,energies,vectors,eigen_workspace);
		gsl_eigen_symmv_sort(energies,vectors,GSL_EIGEN_SORT_VAL_DESC);
		/* H is no more needed : */
		gsl_matrix_free(H);
		gsl_eigen_symmv_free(eigen_workspace);

		tempvector = compute_tempvector(combfactor,vectors,matrix_size);

#ifdef DEBUG
		printf("\nEnergies = \n\n");
		print_vector(energies,matrix_size);
		printf("\ncombfactor = \n\n");
		print_vector_double(combfactor,matrix_size);
		printf("\ntempvector = \n\n");
		print_vector_double(tempvector,matrix_size);
		printf("(inverse)vector = \n\n");
		print_matrix(vectors,matrix_size);
		printf("szmatelem = \n\n");
		print_vector_double(szmatelem,matrix_size);
#endif
		/* Compute points for plotting */
		points = compute_overlap(szmatelem,tempvector,vectors,energies,matrix_size,Np,notpoints,tmax);

		write_to_file(points,Np,M,notpoints);
		/* Free memory */
		free(points);
		free(tempvector);
		gsl_matrix_free(vectors);
		gsl_vector_free(energies);

#ifdef BENCHMARK
		printf("Np=%d completed at : ",Np);
		fflush(stdout);
		system("date");
#endif
		/* Compare results to Mathematica */
		verification(points,M,Np);
	}

	free(combfactor);
	free(szmatelem);
	free(statens);
	free(statensp);
	printf("Done.\n");
	return 0;
}

