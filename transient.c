#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "csparse.h"
#include "transient.h"
#include "options.h"
#include "components.h"
#include "algebra.h"
#include "utility.h"
#include "mna.h"
#include "plot.h"

extern int unique_hash; // this is how many nodes we got hash_table.c

int do_transient = 0;

double tran_step;
double tran_finish;
extern double *G, *C, *dc, *rhs, *TempMatrix, *L, *U; // mna.c
extern int *P;//mna.c
extern double *m;
extern cs *C_s, *G_s;


void transient_analysis_tr()
{
	int size = voltages + inductors + unique_hash;
	int i,j;
	double t, temp, temp2;

	double *swap;
	double *e = (double*) malloc(sizeof(double)*size);
	double *sol  = (double*) malloc(sizeof(double)*size);
	double *rhs0 = (double*) malloc(sizeof(double)*size);

	t = 2/tran_step;

	for (i=0; i<size; i++) {
		for (j=0; j<size; j++ ) {
			c_write(i,j, c_read(i,j)*t); // C[i*size+j] *= t;

			temp = g_read(i,j);
			temp2 = c_read(i,j);

			g_add(i, j, temp2);
			c_write(i, j, temp-temp2);
		}

		sol[i] = dc[i];
		rhs0[i] = rhs[i];
	}

	if ( method_choice == NonIterative ) {
		decompose(size, &P, method_noniter);
	} else {
		if ( sparse_use == 0 ) {
			for (i=0; i<size; i++) {
				m[i] = G[i*size+i];
				if ( fabs(m[i]) < 0.000001 )
					m[i] = 1;
				else
					m[i] = 1/m[i];
			}
		} else {
			cs_get_diag(G_s, m, size);
		}
	}

	if ( sparse_use == 1 )
		G_s = cs_compress(G_s);

	for ( t=tran_step; t <= tran_finish; t+=tran_step ) {
		generate_rhs(rhs, size, unique_hash, 1, t);

		/*for (i=0; i<size; i++) {
			rhs0[i] += rhs[i];
			for (j=0; j<size; j++ )
				rhs0[i] -=  c_read(i,j) * sol[j];
		}*/
		
		for ( i=0; i<size; i++ ) {
			e[i] = 0;
			for ( j=0; j< size; j++ )
				e[i] += c_read(i,j)*sol[j];
			
			e[i]= rhs[i] + rhs0[i] - e[i];
		}

		// rhs = e-1(t) + e(t) -*((G-2/hC)*X(n-1))

		solve(m, P, sol, rhs0, size);
		print_plots(t,sol, P);

		swap = rhs;
		rhs = rhs0;
		rhs0 = swap;
	}

	plot_finalize();

	free(sol);
	free(rhs0);
}

int Doolittle_LU_Decomposition_with_Pivoting(double *A, int pivot[], int n);
int Doolittle_LU_with_Pivoting_Solve(double *A, double B[], int pivot[],
		double x[], int n);


void transient_analysis_be()
{
	int size = voltages + inductors + unique_hash;
	int i,j;
	double t, k;

	double *sol = (double*) malloc(sizeof(double)*size);
	assert(sol);

	t = 1/tran_step;

	for (i=0; i<size; i++) {
		for (j=0; j<size; j++ ) {
			k = c_read(i,j)*t;
			c_write(i, j, k ); // C[i*size+j] *= t;
			g_add(i,j, k); // G[i*size+j] += C[i*size+j];
		}

		sol[i] = 0;
	}


	if ( method_choice == NonIterative ) {
		decompose(size, &P, method_noniter);
	} else {
		if ( sparse_use == 0 ) {
			for (i=0; i<size; i++) {
				m[i] = G[i*size+i];
				if ( fabs(m[i]) < 0.000001 )
					m[i] = 1;
				else
					m[i] = 1/m[i];
			}
		} else {
			cs_get_diag(G_s, m, size);
		}
	}

	if ( sparse_use == 1 ) {
		G_s = cs_compress(G_s);
		assert(G_s);
	}
		

	for ( t=0; t <= tran_finish; t+=tran_step ) {
		generate_rhs(rhs, size, unique_hash, 1, t);

		for (i=0; i<size; i++)
			for (j=0; j<size; j++ )
				rhs[i]+= c_read(i, j) * sol[j];

		// rhs = e(t)+1/h*(C*X(n-1))

		solve(m, P, sol, rhs, size);

		print_array(rhs, size, stdout);
		printf("----\n");
		print_array(sol, size, stdout);
		printf("\n\n");
		print_plots(t,sol, P);
	}

	plot_finalize();

}

void transient_analysis()
{
	if ( method_tran == Tr )
		transient_analysis_tr();
	else if ( method_tran == Be )
		transient_analysis_be();

	printf("[+] Transient analysis: Done\n");
}

