#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "syntactic.tab.h"
#include "solution.h"
#include "utility.h"
#include "csparse.h"

static double *dc_point;
extern int spd_flag;
extern int iter_type;
extern double itol;
extern int use_sparse;
extern struct components_t *g_components;
extern struct instruction_t *g_instructions;
extern struct option_t *g_options;
extern int transient_method;

void instruction_sparse_tran_tr(struct instruction_t *instr, int max_nodes, int sources, 
    cs *MNA_G, cs *MNA_C, double *RHS, double *m)
{
  struct instruction_t *ptr;
  double t;
  int i,j;
  double *solution_0 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *solution_1 = (double*) z_calloc((max_nodes+sources), sizeof(double));
  double *RHS_0 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *RHS_1 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *temp = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  cs *matrix, *matrix2;
  csn *N=NULL;
  css *S=NULL;
  double *swap;

  memcpy(RHS_0, RHS, sizeof(double)*(max_nodes+sources));
  memcpy(solution_0, dc_point, sizeof(double)*(max_nodes+sources));
  
  matrix = cs_add(MNA_G, MNA_C,1.0f,  2/instr->tran.time_step);
  matrix2 = cs_add(MNA_G, MNA_C, 1.0f, -2/instr->tran.time_step);

  if ( iter_type == NoIter ) {
    if ( spd_flag == 0 ) {
      // edw exw LU
      S = cs_sqr (2, matrix, 0) ;              /* ordering and symbolic analysis */
      N = cs_lu (matrix, S, 1) ;                 /* numeric LU factorization */
    } else {
      // edw exw cholesky
      S = cs_schol(1,matrix);
      N = cs_chol(matrix,S);
    }

  } else {
    for (i=0; i<max_nodes+sources; i++ ) 
      m[i] = cs_atxy(matrix, i, i );
  }

  
  for ( t=instr->tran.time_step; t<=instr->tran.time_finish; t+= instr->tran.time_step ) {
    calculate_RHS(g_components, max_nodes, sources, RHS_1, t, 0);
    
    memset(temp,0, sizeof(double)*(max_nodes+sources));
    cs_gaxpy(matrix2, solution_0, temp);
    add_vectors(RHS_0, RHS_1, RHS_0, max_nodes+sources);
    sub_vectors(RHS_0, temp, RHS_0, max_nodes+sources);


    if ( iter_type == NoIter ) {
      if ( spd_flag == 0 ) 
        cs_lusol(S, N, RHS_0, solution_1, (max_nodes+sources));
      else
        cs_cholsol(S,N, RHS_0, solution_1, (max_nodes+sources));
    } else {
      memset(solution_1, 0, sizeof(double) * ( max_nodes+sources));
      if ( iter_type == CG ) {
        conjugate_sparse(matrix, solution_1 , RHS_0, m, itol, max_nodes+sources);
      } else if ( iter_type == BiCG ) {
        biconjugate_sparse(matrix, solution_1 , RHS_0, m, itol, max_nodes+sources);
      }
    }
    
    /*printf("RHS: \n");

    for (i=0; i< max_nodes + sources ;i++ )	
      printf("%7g\n", RHS_1[i]);*/

    /*printf("Result:\n");
    print_array(solution_1,max_nodes+sources);*/

    for ( ptr = g_instructions; ptr!=NULL; ptr=ptr->next) {

      if ( ptr->type == Plot ) {
        for (j=0; j<ptr->plot.num; j++ )
          fprintf( ptr->plot.output[j], "%7G %g\n", 
              t, solution_1[ptr->plot.list[j]]);

      }
    }

    swap = RHS_0;
    RHS_0 = RHS_1;
    RHS_1 = swap;

    swap = solution_0;
    solution_0 = solution_1;
    solution_1 = swap;
  }
  
  cs_spfree(matrix);
  cs_spfree(matrix2);
  free(solution_0);
  free(solution_1);
  free(RHS_0);
  free(RHS_1);
}


void instruction_sparse_tran_be(struct instruction_t *instr, int max_nodes, int sources, 
    cs *MNA_G, cs *MNA_C, double *RHS, double *m)
{
 struct instruction_t *ptr;
  double t;
  int i,j;
  double *solution_0 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *solution_1 = (double*) z_calloc(max_nodes+sources,sizeof(double));
 
  double *RHS_0 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *RHS_1 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *temp = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  cs *matrix ;
  csn *N=NULL;
  css *S=NULL;
  double *swap;


  memcpy(RHS_0, RHS, sizeof(double)*(max_nodes+sources));
  memcpy(solution_0, dc_point, sizeof(double)*(max_nodes+sources));
  
  matrix = cs_add(MNA_G, MNA_C,1.0f,  1/instr->tran.time_step);

  if ( matrix == NULL )
    exit(0);
  
  for (i=0; i<MNA_C->p[ MNA_C->n]; i++ ) {
    MNA_C->x[i] /= instr->tran.time_step;
  }


  if ( iter_type == NoIter ) {
    if ( spd_flag == 0 ) {
      // edw exw LU
      S = cs_sqr (2, matrix, 0) ;              /* ordering and symbolic analysis */
      N = cs_lu (matrix, S, 1) ;                 /* numeric LU factorization */
    } else {
      // edw exw cholesky
      S = cs_schol(1,matrix);
      N = cs_chol(matrix,S);
    }

  } else {
    for (i=0; i<max_nodes+sources; i++ ) 
      m[i] = cs_atxy(matrix, i, i );
  }

  
  for ( t=instr->tran.time_step; t<=instr->tran.time_finish; t+= instr->tran.time_step ) {
    calculate_RHS(g_components, max_nodes, sources, RHS_1, t,0);
    
    memset(temp,0, sizeof(double)*(max_nodes+sources));
    cs_gaxpy(MNA_C, solution_0, temp);
    add_vectors(RHS_1, temp, RHS_0, max_nodes+sources);

    if ( iter_type == NoIter ) {
      if ( spd_flag == 0 ) 
        cs_lusol(S, N, RHS_0, solution_1, (max_nodes+sources));
      else
        cs_cholsol(S,N, RHS_0, solution_1, (max_nodes+sources));
    } else {
      memset(solution_1, 0, sizeof(double) * ( max_nodes+sources));
      if ( iter_type == CG ) {
        conjugate_sparse(matrix, solution_1 , RHS_0, m, itol, max_nodes+sources);
      } else if ( iter_type == BiCG ) {
        biconjugate_sparse(matrix, solution_1 , RHS_0, m, itol, max_nodes+sources);
      }
    }
    
    /*printf("RHS: \n");

    for (i=0; i< max_nodes + sources ;i++ )	
      printf("%7g\n", RHS_0[i]);

    printf("Result:\n");
    print_array(solution_1,max_nodes+sources); */

    for ( ptr = g_instructions; ptr!=NULL; ptr=ptr->next) {

      if ( ptr->type == Plot ) {
        for (j=0; j<ptr->plot.num; j++ )
          fprintf( ptr->plot.output[j], "%7G %g\n", 
              t, solution_1[ptr->plot.list[j]]);

      }
    }

    swap = RHS_0;
    RHS_0 = RHS_1;
    RHS_1 = swap;

    swap = solution_0;
    solution_0 = solution_1;
    solution_1 = swap;
  }
  
  cs_spfree(matrix);
  free(solution_0);
  free(solution_1);
  free(RHS_0);
  free(RHS_1);
}

void instruction_sparse_tran(struct instruction_t *instr, int max_nodes, int sources, 
    cs *MNA_G, cs *MNA_C, double *RHS, double *m)
{
  if ( transient_method == TR ) {
    instruction_sparse_tran_tr(instr, max_nodes, sources, MNA_G, MNA_C, RHS, m);
  } else if ( transient_method == BE ) {
    instruction_sparse_tran_be(instr, max_nodes, sources, MNA_G, MNA_C, RHS, m);
  } else {
    printf("Unknown transient method :%d\n", transient_method);
  }
}

void circuit_mna_sparse(struct components_t *circuit, cs** MNA_G, cs** MNA_C, int *max_nodes, int *sources,
    int element_types[], int **renamed_nodes)
{
  int transient_analysis = 0;
  struct components_t *s ;
  struct instruction_t *w;
  int max_v_id;
  int  elements, non_zero = 0;

  int inductors = 0;
  // arxika metatrepw tous puknwtes se anoixtokuklwma ( tous afairw apo to kuklwma )
  // kai ta phnia se vraxukuklwma, ta antika8istw me phges tashs 

  s = g_components;

  max_v_id = 1;
  elements = *sources = 0;

  while (s) {
    switch(s->data.type) {
      case V:
        max_v_id = ( s->data.t1.id > max_v_id ? s->data.t1.id : max_v_id );
        if ( s->data.t1.is_ground == 0 )
          (*sources)++;
        non_zero ++;
        break;
      case R:
        if ( s->data.t1.plus == *max_nodes || s->data.t1.minus == *max_nodes)
          non_zero++;
        else
          non_zero+=4;
        elements++;
        break;
      case L:
        inductors++;
        break;
    }

    s = s->next;
  }

  *sources += inductors;
  circuit_rename_nodes(g_components, g_instructions, element_types, renamed_nodes, max_nodes );
  circuit_print(g_components);
  printf("inductors         : %d\n", inductors);
  printf("sources           : %d\n"
         "nodes             : %d\n", *sources-inductors, *max_nodes);
  printf("artificial sources: %d\n", *sources);
  printf("non zero elementrs: %d\n", non_zero);

  for ( w = g_instructions; w; w=w->next ){
    if ( w->type == Tran ) {
      transient_analysis = 1;
      break;
    }
  }

  *MNA_G = cs_spalloc(*max_nodes+*sources, *max_nodes+ *sources, non_zero, 1,1);

  if ( transient_analysis )
    *MNA_C= cs_spalloc(*max_nodes+*sources, *max_nodes+ *sources, non_zero, 1,1);
  // meta ftiaxnoume to panw aristera elements x elements pou einai ta pa8htika stoixeia
  for (s=circuit; s!=NULL ;s=s->next) {
    if ( s->data.type == R ) {
      if ( s->data.t1.plus < *max_nodes) {
        cs_add_to_entry(*MNA_G, s->data.t1.plus, s->data.t1.plus,1/s->data.t1.val);
      } if ( s->data.t1.minus < *max_nodes ) {
        cs_add_to_entry(*MNA_G, s->data.t1.minus, s->data.t1.minus, 1/s->data.t1.val);
      } if ( s->data.t1.minus < *max_nodes && s->data.t1.plus < *max_nodes ) {
        cs_add_to_entry(*MNA_G, s->data.t1.minus, s->data.t1.plus, -1/s->data.t1.val);
        cs_add_to_entry(*MNA_G, s->data.t1.plus, s->data.t1.minus, -1/s->data.t1.val);
      } 
    } else if ( transient_analysis == 1 && s->data.type == C ) {
      if ( s->data.t1.plus < *max_nodes) {
        cs_add_to_entry(*MNA_C, s->data.t1.plus, s->data.t1.plus, s->data.t1.val);
      } if ( s->data.t1.minus < *max_nodes ) {
        cs_add_to_entry(*MNA_C, s->data.t1.minus, s->data.t1.minus, s->data.t1.val);
      } if ( s->data.t1.minus < *max_nodes && s->data.t1.plus < *max_nodes ) {
        cs_add_to_entry(*MNA_C, s->data.t1.minus, s->data.t1.plus, -s->data.t1.val);
        cs_add_to_entry(*MNA_C, s->data.t1.plus, s->data.t1.minus, -s->data.t1.val);
      } 
    } else if ( transient_analysis && s->data.type == L ) {
      cs_entry(*MNA_C, *max_nodes+*sources-inductors +s->data.t1.id, 
          *max_nodes+*sources-inductors+s->data.t1.id, - s->data.t1.val);
    }

    if ( (s->data.type == V && s->data.t1.is_ground == 0) || s->data.type == L ) {
      if ( s->data.t1.plus < *max_nodes ) {
        cs_entry(*MNA_G, *max_nodes + s->data.t1.id, s->data.t1.plus, 1.0 );
        cs_entry(*MNA_G, s->data.t1.plus, *max_nodes + s->data.t1.id, 1.0 );
      }
      if ( s->data.t1.minus < *max_nodes ) {
        cs_entry(*MNA_G, *max_nodes + s->data.t1.id,s->data.t1.minus, -1.0 );
        cs_entry(*MNA_G, s->data.t1.minus, *max_nodes + s->data.t1.id, -1.0 );
      }
    }

  }
  cs_print_formated(*MNA_G, "mna_analysis", *max_nodes+*sources);
  if (transient_analysis)
    cs_print_formated(*MNA_C, "mna_analysis_transient", *max_nodes+*sources);
}



int instruction_dc_sparse(struct instruction_t *instr, int max_nodes, int sources, int renamed_nodes[], 
    double *RHS, css *S, csn *N, cs* MNA, double *m, double *result)
{

  double begin;
  double end;
  double step;
  double dummy;
  int i;
  int original_source_val;
  ;
  struct components_t *s;
  struct instruction_t *ptr = g_instructions;
  int j;


  begin = instr->dc.begin;
  end = instr->dc.end;
  step = instr->dc.inc;

  if(instr->dc.sourceType == Voltage){
    for (s=g_components;s!=NULL; s=s->next) {
      if( s->data.type == V && s->data.t1.id == instr->dc.source){
        i = 0;
        original_source_val = s->data.t1.val;
        for(dummy = begin ; dummy <= end ;i++, dummy = dummy+step){

          s->data.t1.val = dummy;
          //printf("Solving for Voltage Source value %g\n",dummy);
          calculate_RHS(g_components,max_nodes,sources,RHS, -1, 1);
          if ( iter_type == NoIter ) {
            if ( spd_flag == 0 ) 
              cs_lusol(S, N, RHS, result, (max_nodes+sources));
            else
              cs_cholsol(S,N, RHS, result, (max_nodes+sources));
          }else {
            memset(result, 0, sizeof(double) * ( max_nodes+sources));
            if ( iter_type == CG ) {
              conjugate_sparse(MNA, result , RHS, m, itol, max_nodes+sources);
            } else if ( iter_type == BiCG ) {
              biconjugate_sparse(MNA, result , RHS, m, itol, max_nodes+sources);
            }
          }

          /*printf("Result:\n");
          print_array(result,max_nodes+sources);*/

          for ( ptr = g_instructions; ptr!=NULL; ptr=ptr->next) {

            if ( ptr->type == Plot ) {
              for (j=0; j<ptr->plot.num; j++ )
                fprintf( ptr->plot.output[j], "%G %g\n", 
                    dummy, result[ptr->plot.list[j]]);

            }

          }
        }
        s->data.t1.val = original_source_val;
      } 
    }
  } else if  (instr->dc.sourceType == Current ) {
    for (s=g_components;s!=NULL; s=s->next) {
      if( s->data.type == I && s->data.t1.id == instr->dc.source){
        i = 0;
        original_source_val = s->data.t1.val;
        for(dummy = begin ; dummy <= end ;i++, dummy = dummy+step){

          s->data.t1.val = dummy;
          //printf("Solve for Current Source value %g\n",dummy);

          calculate_RHS(g_components,max_nodes,sources,RHS, -1, 1);
          if ( iter_type == NoIter ) {
            if ( spd_flag == 0 ) 
              cs_lusol(S, N, RHS, result, (max_nodes+sources));
            else
              cs_cholsol(S,N, RHS, result, (max_nodes+sources));
          }else {
            memset(result, 0, sizeof(double) * ( max_nodes+sources));
            if ( iter_type == CG ) {
              conjugate_sparse(MNA, result , RHS, m, itol, max_nodes+sources);
            } else if ( iter_type == BiCG ) {
              biconjugate_sparse(MNA, result , RHS, m, itol, max_nodes+sources);
            }
          }

          /*printf("Result:\n");
          print_array(result,max_nodes+sources);*/

          for (ptr = g_instructions; ptr!=NULL; ptr=ptr->next) {

            if ( ptr->type == Plot ) {
              for (j=0; j<ptr->plot.num; j++ )
                fprintf( ptr->plot.output[j], "%G %g\n", 
                    dummy, result[ptr->plot.list[j]]);

            }

          }}
          s->data.t1.val = original_source_val;
      } 
    }
  }

  /*for ( i=0; i < max_nodes+sources; i++ ) {
    printf("%d -- %d (%g)\n", S->q[i], N->pinv[i], result[i]);
    }*/
  return 0;

}


int execute_instructions_sparse(cs *MNA_sparse_G, cs *MNA_sparse_C, int max_nodes, int sources, int *renamed_nodes, int stoixeia[])
{
  double *RHS = NULL;
  double *m=NULL;
  int i;
  double *result =NULL;
  struct instruction_t *instr;
  // Gia tous sparse pinakes
  cs *MNA_compressed_G = NULL, *MNA_compressed_C = NULL;
  css *S =  NULL;
  csn *N = NULL;

  instr = g_instructions;

  RHS = (double*) z_calloc(max_nodes+sources, sizeof(double));
  result = (double*) z_calloc((max_nodes+sources), sizeof(double));
  dc_point = (double*) z_calloc((max_nodes+sources), sizeof(double));

  MNA_compressed_G = cs_compress(MNA_sparse_G);

  if ( MNA_sparse_C ) {
    MNA_compressed_C = cs_compress(MNA_sparse_C);
  }
  
  
  cs_spfree(MNA_sparse_G);
  calculate_RHS(g_components,max_nodes,sources,RHS,-1, 1);

  if ( iter_type == NoIter ) {
    // Sparse Matrices

    if ( spd_flag == 0 ) {
      // edw exw LU
      S = cs_sqr (2, MNA_compressed_G, 0) ;              /* ordering and symbolic analysis */
      N = cs_lu (MNA_compressed_G, S, 1) ;                 /* numeric LU factorization */
      cs_lusol(S, N, RHS, dc_point, (max_nodes+sources));
      
    } else {
      // edw exw cholesky
      S = cs_schol(1,MNA_compressed_G);
      N = cs_chol(MNA_compressed_G,S);
      cs_cholsol(S, N, RHS, dc_point, (max_nodes+sources));
    }
  } else {
    m = (double*) z_malloc(sizeof(double) * (max_nodes+sources));

    for (i=0; i<max_nodes+sources; i++ ) 
      m[i] = cs_atxy(MNA_compressed_G, i, i );

    if ( biconjugate_sparse(MNA_compressed_G, dc_point, RHS, m, itol, max_nodes+sources) == -1 ) {
      printf("[-] Biconjugate Failed to solve the DC_point\n");
      exit(0);
    }
  }

  printf("Circuit Solution\n");
  print_array(dc_point, max_nodes+sources);

  calculate_RHS(g_components,max_nodes,sources,RHS,-1, 1);

  FILE *dc = fopen("dc_point", "w");
  for (i=0; i<max_nodes+sources; i++)
    fprintf(dc,"%10g\n", dc_point[i]);
  fclose(dc);


  while ( instr ) {
    if(instr->type == Dc) {
      instruction_dc_sparse(instr, max_nodes, sources, renamed_nodes,
          RHS, S, N, MNA_compressed_G,m, result);
    } else if ( instr->type == Tran ) {
      instruction_sparse_tran(instr, max_nodes, sources, MNA_compressed_G, MNA_compressed_C,
        RHS, m);
    }
    instr = instr->next;
  }

  if (dc_point )
    free(dc_point);

  if (  S  )
    cs_sfree(S);
  if ( N ) 
    cs_nfree(N);

  if ( MNA_compressed_G )
    cs_spfree(MNA_compressed_G);

  if ( MNA_compressed_C )
    cs_spfree(MNA_compressed_C);

  if ( m ) free(m);
  if ( result) free(result);
  if ( RHS ) free(RHS);
  if ( renamed_nodes ) free(renamed_nodes);

  return 0;
}




