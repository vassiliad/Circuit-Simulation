#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "options.h"
#include "csparse.h"
#include "components.h"
#include "algebra.h"

extern double *G, *C;
extern cs *G_s, *C_s;
extern css* S;
extern csn* N;

static double *_r;
static double *_p;
static double *_temp;
static double *_Ap ;
static double *_z = NULL;


void multiply_vector_scalar(double *vec, double scalar, double* output, int size)
{
	int i =0;

	for (i=0; i<size; i++ )
		output[i] = vec[i] * scalar;
}


double dot_vectors(double *v1, double *v2, int size ) {
	int i =0;
	double output = 0;

	for (i=0; i<size; i++ )
		output+= v1[i] * v2[i];
	return output;
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


int Doolittle_LU_Decomposition_with_Pivoting(double *A, int pivot[], int n)
{
	int i, j, k, new_p;
	double *p_k, *p_row, *p_col;
	long double max;

	for ( i=0; i<n; i ++ )
		pivot[i] = i;

	for (k = 0, p_k = A; k < n; p_k += n, k++) {

		new_p = pivot[k];
		max = fabs( *(p_k + k) );
		for (j = k + 1, p_row = p_k + n; j < n; j++, p_row += n) {
			if ( max < fabs(*(p_row + k)) ) {
				max = fabs(*(p_row + k));
				new_p = j;
				p_col = p_row;
			}
		}

		if (pivot[k] != new_p) {
			for (j = 0; j < n; j++) {
				max = *(p_k + j);
				*(p_k + j) = *(p_col + j);
				*(p_col + j) = max;
			}

			max = pivot[new_p];
			pivot[new_p] = pivot[k];
			pivot[k] = max;
		}


		if ( *(p_k + k) == 0.0 ) {
			printf("[-] Failed to decompose LU\n");
			return -1;
		}


		for (i = k+1, p_row = p_k + n; i < n; p_row += n, i++) {
			*(p_row + k) /= *(p_k + k);
		}  


		for (i = k+1, p_row = p_k + n; i < n; p_row += n, i++)
			for (j = k+1; j < n; j++)
				*(p_row + j) -= *(p_row + k) * *(p_k + j);

	}

	return 0;
}

int Choleski_LU_Decomposition(double *A, int n)
{
	int i, k, p;
	double *p_Lk0;                   // pointer to L[k][0]
	double *p_Lkp;                   // pointer to L[k][p]  
	double *p_Lkk;                   // pointer to diagonal element on row k.
	double *p_Li0;                   // pointer to L[i][0]
	double reciprocal;

	for (k = 0, p_Lk0 = A; k < n; p_Lk0 += n, k++) {

		//            Update pointer to row k diagonal element.   

		p_Lkk = p_Lk0 + k;

		//            Calculate the difference of the diagonal element in row k
		//            from the sum of squares of elements row k from column 0 to 
		//            column k-1.

		for (p = 0, p_Lkp = p_Lk0; p < k; p_Lkp += 1,  p++)
			*p_Lkk -= *p_Lkp * *p_Lkp;

		//            If diagonal element is not positive, return the error code,
		//            the matrix is not positive definite symmetric.

		if ( *p_Lkk <= 0.0 ) *p_Lkk = 1.0;

		//            Otherwise take the square root of the diagonal element.

		*p_Lkk = sqrt( *p_Lkk );
		reciprocal = 1.0 / *p_Lkk;

		//            For rows i = k+1 to n-1, column k, calculate the difference
		//            between the i,k th element and the inner product of the first
		//            k-1 columns of row i and row k, then divide the difference by
		//            the diagonal element in row k.
		//            Store the transposed element in the upper triangular matrix.

		p_Li0 = p_Lk0 + n;
		for (i = k + 1; i < n; p_Li0 += n, i++) {
			for (p = 0; p < k; p++)
				*(p_Li0 + k) -= *(p_Li0 + p) * *(p_Lk0 + p);
			*(p_Li0 + k) *= reciprocal;
			*(p_Lk0 + i) = *(p_Li0 + k);
		}  
	}
	return 0;
}
int decompose(int size, int **P, enum NonIterativeMethods type)
{
	int i = 0;

	if ( type == CholDecomp ) {
		if ( voltages || inductors )
			return -1;
	}

	if ( sparse_use == 0 ) {
		switch ( type ) {
			case LUDecomp:
				if ( *P == NULL )
					*P = (int*) malloc(size*sizeof(int));

				return Doolittle_LU_Decomposition_with_Pivoting(G, *P, size);
				break;

			case CholDecomp:
				if ( *P == NULL )
					*P = (int*) malloc(size* sizeof(int));

				for ( i=0; i<size; i++ )
					(*P)[i] = i;

				return Choleski_LU_Decomposition(G, size);
				break;

			default:
				return -2;
		}
	} else {
		cs *comp;

		comp = cs_compress(G_s);

		assert(comp);

		switch( type ) {
			case LUDecomp:
				S = cs_sqr (2, comp, 0) ;              /* ordering and symbolic analysis */
				N = cs_lu (comp, S, 1) ;                 /* numeric LU factorization */
				break;

			case CholDecomp:
				S = cs_schol(1,comp);
				N = cs_chol(comp,S);
				break;
		}

		if ( !S || !N )
			return -1;
	}
	return 0;
}

int biconjugate_sparse(cs *A, double *x, double *b, double *m, double itol, int size){
	double rsold;
	double rsnew;
	double alpha;
	double beta;

	int i;

	memset(_z, sizeof(double)*size, 0);
	memset(_r, sizeof(double)*size, 0);
	memset(_temp, sizeof(double)*size, 0);
	memset(_Ap, sizeof(double)*size, 0);
	memset(_p, sizeof(double)*size, 0);

	//multiply_matrix_vector(A, x, _p, size);
	cs_gaxpy(A, x, _p); // gia na doulepsei swsta prepei _p=[0,...,0]
	sub_vectors(b, _p, _r, size);

	for( i = 0 ; i < size; i++){
		multiply_vector_vector(_r,m,_z,size);
		rsold = rsnew;
		rsnew = dot_vectors(_z,_r,size);
		if(!rsnew){
			printf("biconjugate FAILS \n");
			return -1;
		}
		if(i == 0 ){
			memcpy(_p, _z, size*sizeof(double));
		}
		else{
			beta = rsnew/rsold;
			multiply_vector_scalar(_p,beta,_p,size);
			add_vectors(_z,_p,_p,size);
		}
		
		memset(_Ap, sizeof(double)*size, 0);
		cs_gaxpy(A, _p, _Ap);
		//multiply_matrix_vector(A,_p,_Ap,size);
		alpha = rsnew/dot_vectors(_p,_Ap,size);

		multiply_vector_scalar(_p,alpha,_temp,size);
		add_vectors(x,_temp,x,size);

		multiply_vector_scalar(_Ap,alpha,_temp,size);
		sub_vectors(_r,_temp,_r,size);

		if ( sqrt(rsnew) < itol )
			break;

	}

	return 0;
}

void conjugate_sparse(cs *A, double *x, double *b, double *m, double itol, int size) {
	double rsold;
	double rsnew;
	double alpha;

	int i;

	memset(_z, sizeof(double)*size, 0);
	memset(_r, sizeof(double)*size, 0);
	memset(_temp, sizeof(double)*size, 0);
	memset(_Ap, sizeof(double)*size, 0);
	memset(_p, sizeof(double)*size, 0);
	
	cs_gaxpy(A, x, _p);
	//multiply_matrix_vector(A, x, _p, size);
	sub_vectors(b, _p, _r, size);
	multiply_vector_vector(m, _r, _z, size);

	memcpy(_p, _z, size*sizeof(double));
	rsold = dot_vectors(_r,_z,size);

	for (i=0; i<size; i++ ) {
		memset(_Ap, sizeof(double)*size, 0);
		cs_gaxpy(A, _p, _Ap);
		// multiply_matrix_vector(A,_p, _Ap, size);
		alpha = rsold/dot_vectors(_p, _Ap, size);

		multiply_vector_scalar(_p, alpha, _temp, size);
		add_vectors(x, _temp, x, size);

		multiply_vector_scalar(_Ap, alpha, _temp, size);
		sub_vectors(_r, _temp, _r, size);


		multiply_vector_vector(m,_r,_z,size);
		rsnew = dot_vectors(_z,_r, size);

		if ( sqrt(rsnew) < itol )
			break;

		multiply_vector_scalar(_p, rsnew/rsold, _p, size);
		add_vectors(_z, _p, _p, size);
		rsold = rsnew;
	}
}

int biconjugate(double *A, double *x, double *b, double *m, double itol, int size){
	double rsold;
	double rsnew;
	double alpha;
	double beta;

	int i;

	memset(_z, sizeof(double)*size, 0);
	memset(_r, sizeof(double)*size, 0);
	memset(_temp, sizeof(double)*size, 0);
	memset(_Ap, sizeof(double)*size, 0);
	memset(_p, sizeof(double)*size, 0);

	multiply_matrix_vector(A, x, _p, size);
	sub_vectors(b, _p, _r, size);

	for( i = 0 ; i < size; i++){
		multiply_vector_vector(_r,m,_z,size);
		rsold = rsnew;
		rsnew = dot_vectors(_z,_r,size);
		if(!rsnew){
			printf("biconjugate FAILS \n");
			return -1;
		}
		if(i == 0 ){
			memcpy(_p, _z, size*sizeof(double));
		}
		else{
			beta = rsnew/rsold;
			multiply_vector_scalar(_p,beta,_p,size);
			add_vectors(_z,_p,_p,size);
		}

		multiply_matrix_vector(A,_p,_Ap,size);
		alpha = rsnew/dot_vectors(_p,_Ap,size);

		multiply_vector_scalar(_p,alpha,_temp,size);
		add_vectors(x,_temp,x,size);

		multiply_vector_scalar(_Ap,alpha,_temp,size);
		sub_vectors(_r,_temp,_r,size);

		if ( sqrt(rsnew) < itol )
			break;

	}

	return 0;
}

void conjugate(double *A, double *x, double *b, double *m, double itol, int size) {
	if ( matrix_symmetric(A, size) == 0 ) {
		printf("O pinakas den einai symmetrikos ( conjugate )\n");
		exit(0);
	}

	double rsold;
	double rsnew;
	double alpha;

	int i;
	memset(x, sizeof(double)*size, 0);
	memset(_z, sizeof(double)*size, 0);
	memset(_r, sizeof(double)*size, 0);
	memset(_temp, sizeof(double)*size, 0);
	memset(_Ap, sizeof(double)*size, 0);
	memset(_p, sizeof(double)*size, 0);

	multiply_matrix_vector(A, x, _p, size);
	sub_vectors(b, _p, _r, size);
	multiply_vector_vector(m, _r, _z, size);

	memcpy(_p, _z, size*sizeof(double));
	rsold = dot_vectors(_r,_z,size);
	printf("$$$$$$$$$$$$$$$$$$$$$$ %lf \n",itol);
	for (i=0; i<size*100000; i++ ) {
		multiply_matrix_vector(A,_p, _Ap, size);
		alpha = rsold/dot_vectors(_p, _Ap, size);

		multiply_vector_scalar(_p, alpha, _temp, size);
		add_vectors(x, _temp, x, size);

		multiply_vector_scalar(_Ap, alpha, _temp, size);
		sub_vectors(_r, _temp, _r, size);


		multiply_vector_vector(m,_r,_z,size);
		rsnew = dot_vectors(_z,_r, size);

		if ( sqrt(rsnew) < itol )
			break;

		multiply_vector_scalar(_p, rsnew/rsold, _p, size);
		add_vectors(_z, _p, _p, size);
		rsold = rsnew;
	}
}

void solve_iter(double *b, double *x, double *m, int size, enum IterativeMethods type)
{
	if ( _z == NULL ) {
		_z = (double* ) malloc(size* sizeof(double));
		_temp = (double* ) malloc(size* sizeof(double));
		_Ap = (double* ) malloc(size* sizeof(double));
		_p = (double* ) malloc(size* sizeof(double));
		_r = (double* ) malloc(size* sizeof(double));
	}

	if ( sparse_use == 0 ) {
		if ( type == BiCG ) {
			biconjugate(G, x, b, m, itol, size);
		} else if (type == CG) {
			conjugate(G, x, b, m, itol, size);
		} else
			assert( 0 && "Invalid Iterative Method");
	} else {
		if ( type == BiCG ) {
			biconjugate_sparse(G_s, x, b, m, itol, size);
		} else if (type == CG) {
			conjugate_sparse(G_s, x, b, m, itol, size);
		} else
			assert( 0 && "Invalid Iterative Method");
	}
}

void solve_lu(int *p, double *b, double *x, int size, enum NonIterativeMethods type)
{
	int i, j;

	if ( sparse_use == 0 ) {
		for (i=0; i<size; i++ ) {
			x[i] = b[p[i]];

			for ( j=0; j<i; j++ )
				x[i] -= x[j] * G[i*size+j];

			if ( type == CholDecomp )
				x[i] /= G[i*size+i];
		} 

		for ( i=size-1; i>=0; i-- ) {
			for ( j=i+1; j<size; j++ )
				x[i] -= x[j] * G[i*size+j];

			x[i] /= G[i*size+i];
		}
	} else {
		if ( type == CholDecomp )
			cs_cholsol(S, N, b, x, size);
		else if ( type == LUDecomp )
			cs_lusol(S, N, b, x, size);
		else
			assert(0);
	}

}

