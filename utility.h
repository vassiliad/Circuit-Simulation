#ifndef UTILITY_H
#define UTILITY_H
#include <stdio.h>

void print_matrix(double *m, int size, FILE *file);
void generate_rhs(double *rhs, int size, int nodes, int transient, double t); 
int print_array(double *A , int size, FILE* file);
void settozero(double *vec,int size);
#endif
