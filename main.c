#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "lib_grover_simulation.h"

int main(int argc, char **argv)
{
	float res = inner_product(10,10,10);
	printf("%.20f\n",res);
	/*int i;
	for(i=0;i<=10;i++){ 
		printf("%f ",2.0 << i);
	}
	*/
	return 0;
}

