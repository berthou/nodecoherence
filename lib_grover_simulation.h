#ifndef LIB_GROVER_SIMULATION_CC
#define LIB_GROVER_SIMULATION_CC

#define FACTORIAL(x) 		\
	do {					\
		if(x<=1)			\
		x=1.0;			\
		else				\
		for(i=x-1;i>1;i--){x=x*i;} \
	}while(0)

float
inner_product(float k,float kx,int np)
{
	int i,min,max;
	float result=0,
		  k_fact,
		  npk_fact,
		  kx_fact,
		  npkx_fact,
		  kn_fact,
		  npkxn_fact,
		  n_fact,
		  nkxk_fact;

	FACTORIAL(k_fact);
	FACTORIAL(npk_fact);
	FACTORIAL(kx_fact);
	FACTORIAL(npkx_fact);
	FACTORIAL(kn_fact);
	FACTORIAL(npkxn_fact);
	FACTORIAL(n_fact);
	FACTORIAL(nkxk_fact);  

	min = (k-kx > 0) ? k-kx : 0; 
	max = (np-kx > k) ? k : np-k;
	printf("%d %d\n",min,max);
	for(i=min;i<=max;i++){ 
		result+=pow(-1,i)/(kn_fact*npkxn_fact*n_fact*nkxk_fact);
	}
	result*=sqrt((k_fact*npk_fact*kx_fact*npkx_fact)/pow(2,np));
return result;
}



#endif
