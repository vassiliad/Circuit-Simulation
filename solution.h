#ifndef SOLUTION_H
#define SOLUTION_H


void biconjugate(double *A, double *x, double *b, double *m, double itol, int size);
void conjugate(double *A, double *x, double *b, double *m, double itol, int size);

int forward_substitution(double *L, double *RHS, double *y, const int *P, int size);
int backward_substitution(double *U, double *y, double *x, int size);

int Choleski_LDU_Decomposition(double *A,double *L, int size);
int LU_decomposition(double *A, double *L , double *U,int *P, int size);
int calculate_transpose(double *input,double *output,int size);
int print_matrix(double *A , int size);
int print_array(double *A , int size);



#endif
