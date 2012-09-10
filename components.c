#include <stdio.h>
#include <stdlib.h>
#include "components.h"
int voltages = 0;
int currents = 0;
int resistors = 0;
int capacitors = 0;
int inductors = 0;

i_t *p_i = NULL;
v_t *p_v = NULL;
r_t *p_r = NULL;
l_t *p_l = NULL;
c_t *p_c = NULL;

void new_v(int plus, int minus, double value, transient_t *transient)
{
  v_t *v;
  v = (v_t*) calloc(1,sizeof(v_t));
  v->id = voltages++;
  v->plus = plus;
  v->minus = minus;
  v->transient = transient;
  v->val = value;

  v->next = p_v;
  p_v = v;

#ifdef VERBOSE
  printf("V%d %d %d %g\n", voltages-1, plus, minus , value);
#endif
}

void new_i(int plus, int minus, double value, transient_t *transient)
{
  i_t *i;
  i = (i_t*) calloc(1,sizeof(i_t));
  i->id = currents++;
  i->plus = plus;
  i->minus = minus;
  i->transient = transient;
  i->val = value;

  i->next = p_i;
  p_i = i;
#ifdef VERBOSE
  printf("I%d %d %d %g\n", currents-1, plus, minus , value);
#endif
}

void new_r(int plus, int minus, double value)
{
  r_t *r;
  r = (r_t*) calloc(1,sizeof(r_t));
  r->id = resistors++;
  r->plus = plus;
  r->minus = minus;
  r->val = value;

  r->next = p_r;
  p_r = r;
#ifdef VERBOSE
  printf("R%d %d %d %g\n", resistors-1, plus, minus , value);
#endif
}

void new_c(int plus, int minus, double value)
{
  c_t *c;
  c = (c_t*) calloc(1,sizeof(c_t));
  c->id = capacitors++;
  c->plus = plus;
  c->minus = minus;
  c->val = value;

  c->next = p_c;
  p_c = c;
#ifdef VERBOSE 
  printf("C%d %d %d %g\n", capacitors-1, plus, minus, value );
#endif
}

void new_l(int plus, int minus, double value)
{
  l_t *l;
  l = (l_t*) calloc(1,sizeof(i_t));
  l->id = inductors++;
  l->plus = plus;
  l->minus = minus;
  l->val = value;

  l->next = p_l;
  p_l = l;
#ifdef VERBOSE
  printf("L%d %d %d %g\n", inductors-1, plus, minus , value);
#endif
}

void cleanup_t1(v_t *root)
{
  v_t *p, *next;

  for ( p=root; p; p=next ) {
    next = p->next;
    if ( p->transient ) {
			if ( p->transient->type == Pwl )
				free(p->transient->tpwl.pairs);

      free(p->transient);
		}
    free(p);
  }
}

void cleanup_t2(r_t *root)
{
  r_t *p, *next;

  for ( p=root; p; p=next ) {
    next = p->next;
    free(p);
  }
}

void components_cleanup()
{
  cleanup_t1(p_v);
  cleanup_t1(p_i);
  cleanup_t2(p_r);
  cleanup_t2(p_c);
  cleanup_t2(p_l);

  inductors = resistors = capacitors = voltages = currents = 0;
}

