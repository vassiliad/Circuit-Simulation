#ifndef ANALYSIS_H

#define ANALYSIS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "syntactic.tab.h"
#include "solution.h"
#include "csparse.h"


void circuit_mna(struct components_t *circuit, double **MNA_G, double **MNA_C, int *max_nodes, int *sources,
    int element_types[], int **renamed_nodes);

void solve(double *L, double *U, double *temp, double *result,
    double *RHS,
    int *P, int max_nodes, int sources, double t);

int instruction_dc(struct instruction_t *instr, int max_nodes, int sources, int renamed_nodes[], 
    double *MNA, double *RHS, double *L, double *U, double *m, int *P,
    double *temp, double *result);

int execute_instructions(double *MNA_G, double *MNA_C,  int max_nodes, int sources, int *renamed_nodes, int stoixeia[]);

#endif