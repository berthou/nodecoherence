#ifndef LIB_GROVER_SIMULATION_H
#define LIB_GROVER_SIMULATION_H

#include "reference_data.h"
#include "inner20.h"

#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

double  *get_innerproduct_pointer(int Np);
int verification(double *points,int M,int Np);
void write_to_file(double *data,int Np,int M,int size);
double overlap(double *szmatelem,double *tempvector,gsl_matrix *inversevectors, gsl_vector *energies, double t,int matrix_size);
double *compute_overlap(double *szmatelem,double *tempvector,gsl_matrix *inversevectors, gsl_vector *energies,int matrix_size,int Np,int notpoints,double tmax);
double *compute_tempvector(double *combfactor, gsl_matrix *inversevector,int matrix_size);
void print_vector_double(double *v,int size);
void print_vector(gsl_vector *v,int size);
void print_matrix(gsl_matrix *H,int size);
double binomial_coeff(int n, int k);
void compute_state_list(int number,int base, int M, int *data);
int compute_matrix_size(int n,int p);
double inner_product(double k,double kx,double np);

#include "lib_grover_simulation.cc"

#endif
