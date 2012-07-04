#ifndef LIB_GROVER_SIMULATION_CC
#define LIB_GROVER_SIMULATION_CC

void print_vector(gsl_vector *v,size_t size)
{
	int i;
	for(i=0;i<size;i++)
	{
			printf("%g\t",gsl_vector_get(v,i));
	}
	printf("\n\n");
}

void print_matrix(gsl_matrix *H,size_t size)
{
	int i,j;
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			printf("%g ",gsl_matrix_get(H, i, j));
		}
		printf("\n");
	}
	printf("\n\n");
}


#define COMPUTE_COMBFACTOR(Np,statens,M,var)											\
	do {																				\
		var = 1/(sqrt(pow(2,Np*M)));													\
		combfactor_product=1;													\
		for(combfactor_iterator=0;combfactor_iterator<M;combfactor_iterator++) \
		{ 																		\
			BINOMIAL_COEFF(Np,statens[combfactor_iterator],combfactor_bin);			\
			combfactor_product*=combfactor_bin;										\
		}																		\
		combfactor_product=sqrt(combfactor_product);								\
		var*=combfactor_product;												\
	}while(0)



#define BINOMIAL_COEFF(n,k,var)				\
	do {										\
		var = 1;								\
		if(n-2*k>0)								\
		{										\
			for(binomial_iterator=n;binomial_iterator>=n-k+1;binomial_iterator--)				\
			var=var*binomial_iterator/(n-binomial_iterator+1);		\
		}										\
		else									\
		{										\
			for(binomial_iterator=n;binomial_iterator>=k+1;binomial_iterator--)						\
			var=var*binomial_iterator/(n-binomial_iterator+1);			\
		}											\
	}while(0)

float binomial_coeff(int n, int k)
{
	float result = 1;
	int i;
	if(n-2*k>0)
	{
		for(i=n;i>=n-k+1;i--)
			result=result*i/(n-i+1);
	}
	else
	{
		for(i=n;i>=k+1;i--)
			result=result*i/(n-i+1);
	}
	return result;
}

int *compute_state_list(int number,int base, int array_length)
{
	/* Convert number 'number' in base 'base'
	 * the result is stored in an array of size 'array_length'
	 * with LSB representation
	 */

	// Allocate result array
	int *result = (int*)malloc(array_length*sizeof(int));

	// Start to fill from the end of the array :
	int index=array_length-1;

	while(number != 0)
	{
		result[index] = (number % base);
		number= (number / base);
		index--;
	}
	return result;
}

size_t compute_matrix_size(n,p)
{
	/*
	 * Computes n^p and return the result as a size_t
	 */
	int i;
	size_t result=n;
	for(i=1;i<p;i++) {
		result*=n;
	}
	return result;
}

#define FACTORIAL(x) 		\
	do {					\
		if(x<=1)			\
		x=1.0;			\
		else				\
		for(fact_iterator=x-1;fact_iterator>1;fact_iterator--){x=x*fact_iterator;} \
	}while(0)


double inner_product(double k,double kx,double np)
{
	// fact_iterator is used in the FACTORIAL macro
	int i,
		start,
		stop,
		fact_iterator;

	double result=0,
		  k_fact    = k,
		  npk_fact  = np-k,
		  kx_fact   = kx,
		  npkx_fact = np-kx,
		  kn_fact,
		  npkxn_fact,
		  n_fact,
		  nkxk_fact;

	FACTORIAL(k_fact);
	FACTORIAL(npk_fact);
	FACTORIAL(kx_fact);
	FACTORIAL(npkx_fact);

	/* start = max(k-kx)
	 * stop = min(k,Np-kx)
	 */
	start = (k-kx > 0) ? k-kx : 0; 
	stop = (k <= np-kx) ? k : np-kx;

	// Sum
	for(i=start;i<=stop;i++){ 

		kn_fact=k-i;
		npkxn_fact=np-kx-i;
		n_fact=i;
		nkxk_fact=i+kx-k;

		FACTORIAL(kn_fact);
		FACTORIAL(npkxn_fact);
		FACTORIAL(n_fact);
		FACTORIAL(nkxk_fact);  

		result+=(pow(-1,i)/(kn_fact*npkxn_fact*n_fact*nkxk_fact));
	}
	// Sum * sqrt
	result*=sqrt((k_fact*npk_fact*kx_fact*npkx_fact)/pow(2,np));

	return result;
}

#endif
