#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "utility.h"
#include "components.h"

#define  epsilon 2.71828183
#define  pi      3.14159265


double linear_interpolate(double x, double x0, double y0, double x1, double y1)
{
  return y0 + ((x-x0)*y1 - (x-x0)*y0)/(x1-x0);
}

double calculate_ac(const transient_t *transient, double t)
{
  int i;

  if ( transient == NULL ) 
    return 0;

  switch ( transient->type ) {
    case Exp:
      if ( t < transient->texp.td1 ) {
        return transient->texp.i1;
      } else if ( t < transient->texp.td2 ) {
        return transient->texp.i1 + ( transient->texp.i2 - transient->texp.i1 )*
            ( 1 - pow(epsilon, -(t-transient->texp.td1)/transient->texp.tc1));
      } else {
        return transient->texp.i1 + ( transient->texp.i2 - transient->texp.i1 )*
            ( pow(epsilon, -(t-transient->texp.td2)/transient->texp.tc2) 
            - pow(epsilon, -(t-transient->texp.td1)/transient->texp.tc1));
      }
    break;

    case Sin:
      if ( t < transient->tsin.td ) {
        return transient->tsin.i1 + transient->tsin.ia * sin( 2*pi * transient->tsin.ph / 360);
      } else {
        return transient->tsin.i1 + transient->tsin.ia 
          * sin(2*pi * transient->tsin.fr * (t - transient->tsin.td ) 
          + 2*pi * transient->tsin.ph / 360) 
          * pow(epsilon, -(t-transient->tsin.td)* transient->tsin.df);
      }

    break;
    case Pulse:
      t-= transient->tpulse.per * (int)( t/transient->tpulse.per);

      if ( t < transient->tpulse.td ) {
        return transient->tpulse.i1;
      } else if ( t  < transient->tpulse.td + transient->tpulse.tr ) {
        return  linear_interpolate(t,
                                  transient->tpulse.td,
                                  transient->tpulse.i1,
                                  transient->tpulse.td+transient->tpulse.tr,
                                  transient->tpulse.i2);
      } else if ( t  < transient->tpulse.td + transient->tpulse.tr + 
                       transient->tpulse.pw) {
        return transient->tpulse.i2;
      } else if ( t  < transient->tpulse.td + transient->tpulse.tr +
                       transient->tpulse.pw + transient->tpulse.tf ) {
        return linear_interpolate(t,
                                  transient->tpulse.td+transient->tpulse.tr+transient->tpulse.pw,
                                  transient->tpulse.i2,
                                  transient->tpulse.td+transient->tpulse.tr+transient->tpulse.pw
                                      + transient->tpulse.tf,
                                  transient->tpulse.i1);
      } else {
        return transient->tpulse.i1;
      }
       
    break;

    case Pwl:
      for ( i=0; i<transient->tpwl.size; i++ ) {
        if ( t > transient->tpwl.pairs[i].t ) {
          if( i == transient->tpwl.size-1 ) {
            return transient->tpwl.pairs[i].i;            
          }

          return linear_interpolate(t,
                      transient->tpwl.pairs[i].t,
                      transient->tpwl.pairs[i].i,
                      transient->tpwl.pairs[i+1].t,
                      transient->tpwl.pairs[i+1].i);
        }
      }

      return transient->tpwl.pairs->i;
    break;

    default:
      printf("Invalid calculate_ac type\n");
      exit(0);
  }

}

void generate_rhs(double *rhs, int size, int nodes, int transient, double t)
{
  v_t *v;
  i_t *i;

  memset(rhs, 0, sizeof(double)*size);
  double temp;
  for (i=p_i; i; i=i->next) {
    if ( transient ==1 && i->transient)
      temp = calculate_ac(i->transient, t) + i->val;
		else
			temp = i->val;

    if ( i->plus > 0 )
      rhs[i->plus-1] -= temp;
    if ( i->minus > 0 )
      rhs[i->minus-1] += temp;
  }

  for (v=p_v; v; v=v->next) {
    if ( transient == 1 && v->transient )
      temp = calculate_ac(v->transient, t) + v->val;
		else
			temp = v->val;

    rhs[v->id + nodes] = temp;
  }
}

void print_matrix(double *m, int size, FILE *file)
{
  int i,j;
  for (i=0; i<size; i++ ) {
    for (j=0; j<size; j++ )
      fprintf(file, "%10g\t", m[i*size+j]);
    fprintf(file, "\n");      
  }
}

int print_array(double *A , int size, FILE* file){
	int j;
	for(j = 0 ; j < size  ; j++){
		fprintf(file, "%g\n",A[j]);
	}
  return 0;
}

