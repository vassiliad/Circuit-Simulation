#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "dc_instruction.h"
#include "mna.h"
#include "plot.h"
#include "options.h"
#include "algebra.h"
#include "utility.h"
#include "components.h"

extern int *P;//mna.c
extern double *m;
extern double *dc, *rhs;
extern int mna_size;
extern int unique_hash;

int do_dc_instruction = 0;
double dc_start;
double dc_stop;
double dc_step;
int dc_is_current;
char * dc_id;


void dc_instruction()
{
	double t;

	
	for ( t = dc_start; t<=dc_stop; t+=dc_step) {
		if ( dc_is_current ) {
			i_t *i;
			
			for ( i=p_i; i; i = i->next )
					if ( strcasecmp(i->string_id, dc_id) == 0 ) {
						i->val = t;
						break;
					}
			if ( i == NULL ) {
				printf("[-] Specified current source does not exit\n");
				exit(1);
			}
				
		} else {
			v_t *v;
			
			for ( v=p_v; v; v = v->next )
					if ( strcasecmp(v->string_id, dc_id) == 0 ) {
						v->val = t;
						break;
					}
			
			if ( v == NULL ) {
				printf("[-] Specified voltage source does nto exist\n");
				exit(1);
			}
		}
		
		generate_rhs(rhs, mna_size, unique_hash, 0, 0);
		solve(m, P, dc, rhs, mna_size);
		print_plots(t, dc, P);
	}
	
	plot_finalize();
}
