#include <math.h>                                     // required for fabs()
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "options.h"
#include "mna.h"
#include "components.h"
#include "algebra.h"
#include "utility.h"
#include "csparse.h"

extern int unique_hash; // this is how many nodes we got

double *G=NULL, *C=NULL, *dc, *rhs, *m=NULL;
cs *G_s=NULL, *C_s=NULL;
css *S=NULL;
csn *N=NULL;

int *P=NULL;
int mna_size=0;


int Doolittle_LU_Decomposition(double *A, int n);
int Doolittle_LU_Solve(double *LU, double B[], double x[], int n);

void mna_free()
{
  if ( G )  free(G);
	if ( C )  free(C);
	if ( P )  free(P);

	free(dc);
	free(rhs);
}

 
 
void solve_dc()
{
	int i, ret;
	cs *temp;
	P = NULL;
	dc = (double*) calloc(mna_size, sizeof(double));
	rhs = (double*) calloc(mna_size*mna_size, sizeof(double));
//	FILE *dc_point;

	generate_rhs(rhs, mna_size, unique_hash, 0, 0);

  if ( method_choice == NonIterative ) {
			ret =  decompose( mna_size, &P, method_noniter);
			printf("asdjkgyasjkdfyasldjavskdjhgas  %d\n\n\n", ret);
			if (ret != 0 ) {
				printf("[-] Circuit doesn't have dc point\n");
				exit(0);
			} else {
				solve_lu( P, rhs, dc, mna_size, method_noniter);
			}
		} else {
			m = (double*) calloc(mna_size, sizeof(double));

			assert(m!=NULL && "Not enough memory for m");
			
			if ( sparse_use == 0 ) {
				for (i=0; i<mna_size; i++) {
					m[i] = G[i*mna_size+i];
					if ( fabs(m[i]) < 0.000001 )
						m[i] = 1;
					else
						m[i] = 1/m[i];
				}
			} else {
				cs_get_diag(G_s, m, mna_size);
				temp = G_s;
				G_s = cs_compress(G_s);
			}
			solve_iter( rhs, dc, m, mna_size, method_iter);
			G_s = temp;
		}

	

	
//	dc_point = fopen("dc_point", "w");
	printf("[#] DC point in file \"%s\"\n",name_of_file);
	print_array(dc, mna_size, f);
	fclose(f);
#ifdef VERBOSE
	print_array(dc, mna_size, stdout);
	printf("[$] Rhs\n");
	print_array(rhs, mna_size, stdout);
#endif
}

void c_add(int row, int col, double val)
{
	if ( sparse_use == 0 ) {
		C[row*mna_size+col] += val;
	} else
		cs_add_to_entry(C_s, row, col, val);
}

void g_add(int row, int col, double val)
{
	if ( sparse_use == 0 ) {
		G[row*mna_size+col] += val;
	} else
    cs_add_to_entry(G_s, row, col, val);
}

double g_read(int row, int col)
{
	if ( sparse_use == 0 ) {
		return G[row*mna_size+col];
	} else
		return cs_atxy(G_s, row, col);
}

double c_read(int row, int col)
{
  if ( sparse_use == 0 ) {
    return C[row*mna_size+col];
  } else {
    return cs_atxy(C_s, row, col);
	}
}

void g_write(int row, int col, double val)
{
	if ( sparse_use == 0 ) {
		G[row*mna_size+col] = val;
	} else
    cs_entry(G_s, row, col, val);
}

void c_write(int row, int col, double val)
{
	if ( sparse_use == 0 ) {
		C[row*mna_size+col] = val;
	} else
    cs_entry(G_s, row, col, val);
		
}
void mna_analysis()
{
	FILE *g_file, *c_file;
	r_t *r;
	v_t *v;
	l_t *l;
	c_t *c;

	unique_hash--;

	mna_size = voltages + inductors +unique_hash;
#ifdef VERBOSE
	printf("[$] Voltages : %d\n"
			"[$] Inductors: %d\n"
			"[$] Nodes    : %d\n"
			"[$] Total    : %d\n",
			voltages, inductors, unique_hash, mna_size);
#endif

	if ( sparse_use == 0 ) {
		G = (double*) calloc(mna_size*mna_size,sizeof(double));
		C = (double*) calloc(mna_size*mna_size, sizeof(double));
    
    assert(G);
    assert(C);
	} else {
    G_s = cs_spalloc(mna_size, mna_size, voltages + resistors*4, 1,1);
    C_s = cs_spalloc(mna_size, mna_size, voltages + resistors*4, 1,1);

    assert(G_s);
    assert(C_s);
  }

	for ( r =p_r; r; r=r->next ) {
		if ( r->plus > 0 )
			g_add(r->plus-1, r->plus-1, 1/r->val);

		if ( r->minus > 0 )
			g_add(r->minus-1, r->minus-1, 1/r->val);

		if ( r->plus >0 && r->minus > 0 ) {
			g_add(r->minus-1, r->plus-1, -1/r->val);
			g_add(r->plus-1, r->minus-1, -1/r->val);
		}
	}

	for ( v = p_v; v; v=v->next ) {
		if ( v->plus > 0 ) {
			g_write(v->plus -1,unique_hash+ v->id, 1);
			g_write(unique_hash + v->id, v->plus-1, 1);
		}

		if ( v->minus > 0 ) {
			g_write(v->minus-1, unique_hash+v->id, -1);
			g_write(unique_hash +v->id, v->minus-1, -1);
		}
	}

	for ( l = p_l; l; l=l->next ) {
		if ( l->plus > 0 ) {
			g_write(l->plus -1,unique_hash+ l->id + voltages, 1);
			g_write(unique_hash + l->id+voltages, l->plus-1, 1);
		}

		if ( l->minus > 0 ) {
			g_write(l->minus-1, unique_hash+l->id+voltages, -1);
			g_write(unique_hash +l->id+voltages, l->minus-1, -1);
		}

		c_write(voltages+unique_hash+l->id, voltages+unique_hash+l->id, -l->val);
	}

	for ( c=p_c; c; c=c->next) {
		if ( c->plus > 0 )
			c_add(c->plus-1, c->plus-1, c->val);

		if ( c->minus > 0 )
			c_add(c->minus-1, c->minus-1,c->val);

		if ( c->plus >0 && c->minus > 0 ) {
			c_add(c->minus-1, c->plus-1, c->val);
			c_add(c->plus-1, c->minus-1, c->val);
		}
	}

	g_file = fopen("G_matrix", "w");
	if ( sparse_use == 0 )
		print_matrix(G, mna_size, g_file);
	else
    cs_print_formated(G_s, g_file, mna_size);
	fclose(g_file);

	c_file = fopen("C_matrix", "w");
	if ( sparse_use == 0 )
		print_matrix(C, mna_size, c_file);
	else
    cs_print_formated(C_s, c_file, mna_size);
	fclose(c_file);

	printf("[#] G and C are saved in G_matrix and C_matrix\n");
}

