#ifndef LIB_GROVER_SIMULATION_CC
#define LIB_GROVER_SIMULATION_CC

void write_to_file(double *data,int Np,int M,int size)
{
	int i,
		file_desc;

	char x[20],
		 y[20];

	char filename[15],
		 m[2],
		 n[5];

	/* Build filename = "M(value of M)N(value of N).dat"
	 */
	memset(filename,0,15);
	memset(m,0,2);
	memset(n,0,5);
	strcat(filename,"M");
	snprintf(m,sizeof(m),"%d",M);
	strcat(filename,m);
	strcat(filename,"N");
	snprintf(n,sizeof(n),"%d",Np);
	strcat(filename,n);
	strcat(filename,".dat");
	
	/* Open file */
	if((file_desc=open(filename,O_WRONLY | O_CREAT | O_TRUNC, 00600))== -1)
	{
		printf("Unable to open the file %s\n",filename);
		return;
	}

	/* Write file : */
	for(i=0;i<size+1;i++)
	{
		snprintf(x,sizeof(x),"%.15f\t",data[i]);
		write(file_desc,x,strlen(x));
		snprintf(y,sizeof(y),"%.15f\n",data[i+size+1]);
		write(file_desc,y,strlen(y));
	}
	close(file_desc);
}

double overlap(double *szmatelem,double *tempvector,gsl_matrix *inversevectors, gsl_vector *energies, double t,int matrix_size)
{
	int nstilde,
		ms;
	double var=0.0;
	gsl_complex c_sum;

	for(nstilde=0;nstilde<matrix_size;nstilde++)
	{
		c_sum=gsl_complex_polar(0,0);
		for(ms=0;ms<matrix_size;ms++)
		{
			c_sum=gsl_complex_add(c_sum,gsl_complex_polar(tempvector[ms]*gsl_matrix_get(inversevectors,nstilde,ms),(-1)*t*gsl_vector_get(energies,ms)));
		}
		var+=szmatelem[nstilde]*gsl_complex_abs2(c_sum);
	}
	return var;

}

double *compute_overlap(double *szmatelem,double *tempvector,gsl_matrix *inversevectors, gsl_vector *energies,int matrix_size,int Np,int notpoints,double tmax)
{
	int i;
	double *result = (double*)malloc(2*(notpoints+1)*sizeof(double));
	for(i=0;i<notpoints+1;i++)
	{
		result[i]=i*(tmax/notpoints);
		result[i+notpoints+1]=overlap(szmatelem,tempvector,inversevectors,energies,((double)i*(tmax/notpoints))/Np,matrix_size);
		/* printf("t= %f\tSz/Nt= %.10f\n",result[i],result[i+notpoints+1]); */
	}
	return result;
}



inline double *compute_tempvector(double *combfactor, gsl_matrix *inversevector,int matrix_size)
{
	double *result = (double*)malloc(matrix_size*sizeof(double));
	int i,j;
	for(i=0;i<matrix_size;i++)
	{
		result[i]=0.0;
		for(j=0;j<matrix_size;j++)
		{
			result[i]+=combfactor[j]*gsl_matrix_get(inversevector,j,i);
		}
	}
	return result;
}


void print_vector_double(double *v,int size)
{
	int i;
	for(i=0;i<size;i++)
	{
		printf("%.10f\t",v[i]);
	}
	printf("\n");
}

void print_vector(gsl_vector *v,int size)
{
	int i;
	for(i=0;i<size;i++)
	{
		printf("%0.10g\t",gsl_vector_get(v,i));
	}
	printf("\n");
}

void print_matrix(gsl_matrix *H,int size)
{
	int i,j;
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			printf("%0.10g\t",gsl_matrix_get(H, i, j));
		}
		printf("\n");
	}
	printf("\n");
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

double binomial_coeff(int n, int k)
{
	double result = 1;
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

double *compute_state_list(int number,int base, int array_length)
{
	/* Convert number 'number' in base 'base'
	 * the result is stored in an array of size 'array_length'
	 * with LSB representation
	 */

	/* Allocate result array */
	double *result = (double*)malloc(array_length*sizeof(double));

	/* Start to fill from the end of the array : */
	int index=array_length-1;

	while(number != 0)
	{
		result[index] = (number % base);
		number= (number / base);
		index--;
	}
	return result;
}

int compute_matrix_size(int n,int p)
{
	/*
	 * Computes n^p and return the result as a int
	 */
	int i;
	int result=n;
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
		for(fact_iterator=x-1;fact_iterator>1;fact_iterator--){x=(x*fact_iterator);} \
	}while(0)


double inner_product(double k,double kx,double np)
{
	/* fact_iterator is used in the FACTORIAL macro */
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

	/* Sum */
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
	/* Sum * sqrt */
	result*=sqrt((k_fact*npk_fact*kx_fact*npkx_fact)/pow(2,np));

	return result;
}

#endif
