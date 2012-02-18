#ifndef SPARSE_ANALYSIS_H
#define SPARSE_ANALYSIS_H

void circuit_mna_sparse(struct components_t *circuit, cs** MNA_G, cs** MNA_C, int *max_nodes, int *sources,
    int element_types[], int **renamed_nodes);


int instruction_dc_sparse(struct instruction_t *instr, int max_nodes, int sources, int renamed_nodes[], 
    double *RHS, css *S, csn *N, cs* MNA, double *m, double *result);

int instruction_dc(struct instruction_t *instr, int max_nodes, int sources, int renamed_nodes[], 
    double *MNA, double *RHS, double *L, double *U, double *m, int *P,
    double *temp, double *result);

int execute_instructions_sparse(cs *MNA_sparse_G, cs *MNA_sparse_C, int max_nodes, int sources, int *renamed_nodes, int stoixeia[]);

#endif