#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "syntactic.tab.h"
#include "solution.h"
#include "csparse.h"
#include "utility.h"
#include "sparse_analysis.h"
#include "analysis.h"

int yylex_destroy();
int yyparse();
extern struct components_t *g_components;
extern struct instruction_t *g_instructions;
extern struct option_t *g_options;
int spd_flag = 0;
int iter_type = NoIter;
double itol = 0.001;
int use_sparse = 0;
double *dc_point=NULL;

extern FILE * yyin;

void options_cleanup(struct option_t *g_options);
int circuit_rename_nodes(struct components_t *circuit, struct instruction_t *instr, int element_types[], int **renamed_nodes, int *max_nodes);




int main(int argc, char* argv[])
{
  int ret;
  double *MNA_G=NULL, *MNA_C=NULL;
  cs *MNA_sparse_G=NULL, *MNA_sparse_C = NULL;
  int max_nodes, sources, *renamed_nodes;
  int stoixeia[8] = { 0,0,0,0,0,0,0,0 };
  struct option_t *o;
  use_sparse = 0 ;
  g_instructions = NULL;
  g_components = NULL;
  g_options = NULL;

  if ( argc != 2 ) {
    print_help(argv[0]);
    return 0;
  }

  yyin = fopen(argv[1], "r");

  if ( yyin == NULL ) {
    printf("[-] Could not open input file \"%s\"\n", argv[1]);
    return 1;
  }

  ret = yyparse();
  fclose(yyin);
  if ( ret == 0 )
    printf("[+] No errors\n");
  else
    return 1;

  circuit_print(g_components);
  instructions_print(g_instructions);

  for (o=g_options; o!=NULL; o=o->next) {
    switch ( o->type ) {
      case SPD:
        spd_flag = 1;
        break;

      case ITER:
        iter_type = o->iter_type;
        break;

      case ITOL:
        itol = o->itol;
        break;

      case SPARSE:
        use_sparse = 1;
        break;
    }
  }


  if ( use_sparse==0 ){
    circuit_mna(g_components,&MNA_G, &MNA_C,&max_nodes,&sources, stoixeia, &renamed_nodes);
    execute_instructions(MNA_G, MNA_C, max_nodes,
      sources, renamed_nodes, stoixeia);    
  }
  else{
    circuit_mna_sparse(g_components,&MNA_sparse_G, &MNA_sparse_C,&max_nodes,
      &sources, stoixeia, &renamed_nodes);
      execute_instructions_sparse(MNA_sparse_G, MNA_sparse_C, max_nodes,
      sources, renamed_nodes, stoixeia);

  }
  
  circuit_cleanup(g_components);
  instructions_cleanup(g_instructions);
  options_cleanup(g_options);

  yylex_destroy();
  return 0;
}

