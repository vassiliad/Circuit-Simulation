#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

int transpose(double *input,double *output,int size){
	int i,j;
	for(i = 0 ; i < size ; i++)
		for( j = 0; j < size ; j++)
			output[j*size+i] = input[i*size+j];
  return 0;
}

int matrix_symmetric(double *mat, int size)
{
	int i = 0;
	int j = 0;

	for ( i=0; i<size; i++ ) 
		for ( j=0; j<i; j++ )
			if ( fabs(mat[i*size+j]- mat[j*size+i]) > 1e-10 )
				return 0;
	return 1;
}

void multiply_matrix_matrix(double *a, double *b, double *c, int size)
{
  int i ,j , k;

  for ( i=0; i  < size; i ++ ) {
    for ( j=0; j < size; j ++ ) {
      c[i *size + j ] = 0;
      for ( k =0 ; k < size; k ++ ) {
        c[i*size+j] += a[i*size +k] * b[k*size+j];
      }
    }
  }
}

void multiply_matrix_vector(double *mat, double *vector, double *output, int size )
{
  int i, j;

  for (i=0; i<size; i++ ) {
    output[i] = 0;

    for (j=0; j<size; j++ )
      output[i]+= mat[i*size+j] * vector[j];
  }
}

void multiply_vector_vector(double *v1, double *v2, double *output, int size) 
{
  int i=0;
  
  for (i=0; i<size; i++ )
    output[i] = v1[i] * v2[i];
}

void sub_vectors(double *v1, double *v2, double *output, int size) {
  int i=0;
  
  for (i=0; i<size; i++ )
    output[i] = v1[i] - v2[i];
}

void add_vectors(double *v1, double *v2, double *output, int size) {
  int i=0;
  
  for (i=0; i<size; i++ )
    output[i] = v1[i] + v2[i];
}

double dot_vectors(double *v1, double *v2, int size ) {
  int i =0;
  double output = 0;

  for (i=0; i<size; i++ )
    output+= v1[i] * v2[i];
  return output;
}

void multiply_vector_scalar(double *vec, double scalar, double* output, int size)
{
  int i =0;
  
  for (i=0; i<size; i++ )
    output[i] = vec[i] * scalar;
}

void biconjugate(double *A, double *x, double *b, double *m, double itol, int size){
  double *r = (double*) calloc(size, sizeof(double));
  double *rtilde = (double*) calloc(size, sizeof(double));
  double *p = (double*) calloc(size, sizeof(double));
  double *ptilde = (double*) calloc(size, sizeof(double));
  double *temp = (double*) calloc(size, sizeof(double));
  double rsold;
  double *Ap = (double*) calloc(size, sizeof(double));
  double *Aptilde = (double*) calloc(size, sizeof(double));
  double *z = (double* ) calloc(size, sizeof(double));
  double *ztilde = (double* ) calloc(size, sizeof(double));
  double rsnew;
  double alpha;
  double beta;
  
  double *Atilde = (double* ) calloc(size*size, sizeof(double));
  
  int i;
  
  multiply_matrix_vector(A, x, p, size);
  sub_vectors(b, p, r, size);
  memcpy(rtilde, r, size*sizeof(double));
  
  transpose(A,Atilde,size);
  
  for( i = 0 ; i < size; i++){
    multiply_vector_vector(r,m,z,size);
    multiply_vector_vector(rtilde,m,ztilde,size); //epilush tou M^T*z = r
    rsold = rsnew;
    rsnew = dot_vectors(z,rtilde,size);
    if(!rsnew){
      printf("biconjugate FAILS \n");
      return;
    }
    if(i == 0 ){
      memcpy(p, z, size*sizeof(double));
      memcpy(ptilde, ztilde, size*sizeof(double));
    }
    else{
      beta = rsnew/rsold;
      multiply_vector_scalar(p,beta,p,size);
      add_vectors(z,p,p,size);
      multiply_vector_scalar(ptilde,beta,ptilde,size);
      add_vectors(ztilde,ptilde,ptilde,size);
    }
    
    multiply_matrix_vector(A,p,Ap,size);
    multiply_matrix_vector(Atilde,ptilde,Aptilde,size);
    alpha = rsnew/dot_vectors(ptilde,Ap,size);
    
    multiply_vector_scalar(p,alpha,temp,size);
    add_vectors(x,temp,x,size);
   
    multiply_vector_scalar(Ap,alpha,temp,size);
    sub_vectors(r,temp,r,size);
    
    multiply_vector_scalar(Aptilde,alpha,temp,size);
    sub_vectors(rtilde,temp,rtilde,size);
   
     if ( sqrt(rsnew) < itol )
      break;
    
  }
  
  free(r);
  free(rtilde);
  free(p);
  free(ptilde);
  free(temp);
  free(Ap);
  free(Aptilde);
  free(z);
  free(ztilde);
  free(Atilde);
  
}




void conjugate(double *A, double *x, double *b, double *m, double itol, int size) {
	if ( matrix_symmetric(A, size) == 0 ) {
		printf("O pinakas den einai symmetrikos ( conjugate )\n");
		exit(0);
	}

  double *r = (double*) calloc(size, sizeof(double));
  double *p = (double*) calloc(size, sizeof(double));
  double *temp = (double*) calloc(size, sizeof(double));
  double rsold;
  double *Ap = (double*) calloc(size, sizeof(double));
  double *z = (double* ) calloc(size, sizeof(double));
  double rsnew;
  double alpha;
  
  int i;
  
  multiply_matrix_vector(A, x, p, size);
  sub_vectors(b, p, r, size);
  multiply_vector_vector(m, r, z, size);
  
  memcpy(p, z, size*sizeof(double));
  rsold = dot_vectors(r,z,size);
  
  for (i=0; i<size; i++ ) {
    multiply_matrix_vector(A,p, Ap, size);
    alpha = rsold/dot_vectors(p, Ap, size);
    
    multiply_vector_scalar(p, alpha, temp, size);
    add_vectors(x, temp, x, size);
    
    multiply_vector_scalar(Ap, alpha, temp, size);
    sub_vectors(r, temp, r, size);
    

    multiply_vector_vector(m,r,z,size);
    rsnew = dot_vectors(z,r, size);

    if ( sqrt(rsnew) < itol )
      break;
    
    multiply_vector_scalar(p, rsnew/rsold, p, size);
    add_vectors(z, p, p, size);
    rsold = rsnew;
  }
	

  free(r);
  free(p);
  free(temp);
  free(Ap);
  free(z);
}


int backward_substitution(double *U, double *y, double *x, int size)
{
	int i, j;

	double sum;

	x[size-1] = y[size-1] / U[(size-1)*size + (size-1)];

	for (i=size-2; i>=0; i-- ) {
		sum = 0;
		for (j=size-1; j>i; j-- )
			sum += U[i*size+j] * x[j];
		
		x[i] = (y[i]-sum) / U[i*size+i];
	}

	return 1;
}



int forward_substitution(double *L, double *RHS, double *y, const int *P, int size)
{
	int i=0;
	int j;

	double sum;
	
	y[0] = RHS[ P[0] ] / L[0];

	for (i=1; i<size; i++) {
		sum = 0;
		for (j=0; j<i; j++ )
			sum += L[i*size+j] * y[j];
		y[i] = (RHS[P[i]] - sum)/L[i*size+i] ;
	}
	return 1;
}


int Choleski_LDU_Decomposition(double *A,double *L, int size)
{
	int i,j,k;
	double sum;

	for(i = 0 ; i < size ; i++){
		for(j = 0 ; j < size  ; j++){
			L[i*size+j] = 0.0;
		}
	}

	for(i = 0 ; i < size ; i++){

		for (j=0; j<i; j++ ) {
	
			sum = A[i*size+j];
			for (k=0; k<j; k++)
				sum -= L[i*size+k] * L[j*size+k];
			sum /= L[j*size+j];

			L[i*size+j] = sum;
		}
		
		sum = 0.0;
		for (k=0; k<i; k++ )
			sum -= L[i*size+k] * L[i*size+k];
		sum += A[i*size+i];
		
		if ( sum < 0 )	// upo8etoume oti o elegxos gia cholesky exei ginei eksw						
			sum = 0;		// epeidh gnwrizoume oti an exw mono C,R,I to kuklwma mou einai 8etika orismeno

		L[i*size +i] = sqrt( sum );
	}
	return 1;
}

int LU_decomposition(double *A, double *L , double *U,int *P, int size)
{

	int i,j;
	int c;
	int updated;

	long double sum;
	long double temp;
	
	if ( !(A&&L&&U&&P) ) {
		printf("No Memory allocated\n");
		return 0;
	}

	for (i=0; i<size; i++ ) {
		temp = fabs(A[i*size+i]);
		updated = 0;
		for (j=i+1; j<size; j++) {
			if ( temp < fabs(A[j*size+i])) {
				temp = fabs(A[j*size+i]);
				P[i] = j;
				updated = 1;
			}
		}

		if ( updated ) {
			P[P[i]] = i;
			for (j=0; j<size; j++ ) {
				temp = A[i*size +j];
				A[i*size+j] = A[ P[i]*size+j ];
				A[P[i]*size+j] = temp;
			}
		}
	}

	for (i=0; i<size; i++ ) {
		for ( j=0; j<i; j++ ) {
			sum = 0;

			for (c =0; c<j; c++ )
				sum+= L[i*size+c] * U[c*size+j];

			L[i*size+j] = A[i*size+j] - sum;
			L[i*size+j]/=U[j*size+j];
		}

		L[i*size+i] = 1;

		for (j=0; j<size; j++ ) {
			sum = 0;
			for (c=0; c<i; c++ )
				sum+= L[i*size+c] * U[c*size+j];

			U[i*size+j] = A[i*size+j] - sum;
		}
	}

	return 1;
}


int calculate_transpose(double *input,double *output,int size){
	int i,j;
	for(i = 0 ; i < size ; i++)
		for( j = 0; j < size ; j++)
			output[j*size+i] = input[i*size+j];
  return 0;
}



int print_matrix(double *A , int size)
{
	int i,j;
	for(i = 0 ; i < size ; i++){
		for(j = 0 ; j < size  ; j++){
			printf("%.3G\t",A[i*size+j]);
		}
		printf("\n");
	}
  return 0;
}

int print_array(double *A , int size){
	int j;
	for(j = 0 ; j < size  ; j++){
		printf("%g\n",A[j]);
	}
  return 0;
}
/*
	 int main(int argc ,int argv){
	 double A[] = {25,-5,10,-5,17,10,10,10,62};
	 double K[] = {2,1,1,4,-6,0,-2,7,2};
	 double M[] = {5, 1.2, 0.3, -0.6, 1.2, 6, -0.4, 0.9, 0.3, -0.4, 8, 1.7, -0.6, 0.9, 1.7, 10};
	 double U[9];
	 double L[9];
	 int P[3];
	 double result[3];
	 double solution[] = {5,-2,9};
	 memset(P,-1,9*sizeof(int));
	 memset(U,0,9*sizeof(double));
	 memset(L,0,9*sizeof(double));
	 if(solve_system(K,L,U,P,result,solution,3,0)){
	 print_array(result,3);
	 printf("SUCCESS!\n");
	 }
	 else
	 printf("failure\n");
	 }

 */


