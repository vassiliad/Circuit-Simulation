#ifndef UTILITY_H

#define UTILITY_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "syntactic.tab.h"
#include "solution.h"
#include "csparse.h"


#define  epsilon 2.71828183
#define  pi      3.14159265

enum elementTypes {typeV, typeI, typeR, typeC, typeL, typeM, typeD, typeQ };



double linear_interpolate(double x, double x0, double y0, double x1, double y1);
double calculate_ac(const transient_spec_t *transient, double t);
void options_cleanup(struct option_t *g_options);
void instructions_cleanup(struct instruction_t *instr);
void circuit_print(struct components_t *circuit);
int calculate_RHS(struct components_t *circuit,int max_nodes,int sources, double *RHS, double t);
void print_help(char *path);
void instructions_print(struct instruction_t *instr);
void circuit_cleanup(struct components_t *circuit);
int circuit_rename_nodes(struct components_t *circuit, struct instruction_t *instr, int element_types[], int **renamed_nodes, int *max_nodes);

#endif