#ifndef ALGEBRA_H
#define ALGEBRA_H
#include "options.h"

int  decompose(int size,  int **p, enum NonIterativeMethods type);
void solve(double *m , int *P, double *sol, double *rhs,int  size);
void solve_lu(int *p, double *b, double *x,  int size, enum NonIterativeMethods type);
void solve_iter(double *b, double *x, double *m, int size, enum IterativeMethods type);

#endif
