#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int argc , char *argv[]){
	int i,size1,size2;
	double sum , var;
	if (argc!=5){
		printf("./executable file1 file2 size1 size2\n");
		exit(0);
	}
	FILE *f1,*f2;
	f1 = fopen(argv[1],"r");
	f2 = fopen(argv[2],"r");
	size1 = atoi(argv[3]);
	size2 = atoi(argv[4]);
	printf("SIZE %d----%d\n",size1,	size2);
	if(size1!=size2){
		printf("differnet sizes error\n");
		return(0);
	}
	double *first = (double*) calloc (size1,sizeof(double));
	double *second = (double*) calloc (size1,sizeof(double));
	sum = 0.0;
	var = 0.0;
	for(i = 0 ; i < size1 ; i++){
		fscanf(f1," %lf",&first[i]);
		fscanf(f2," %lf",&second[i]);
		sum = sum + fabs(first[i] - second[i]);
		var = var + (first[i] - second[i])*(first[i] - second[i]);
	}
	sum = sum / ((double) size1);
	var = var/  ((double) size1);
	printf("average = %g variance = %g \n",sum,var);
	
	
}