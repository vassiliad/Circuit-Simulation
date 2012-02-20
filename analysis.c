#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "syntactic.tab.h"
#include "solution.h"
#include "csparse.h"
#include "utility.h"


extern int spd_flag;
extern int iter_type;
extern double itol;
extern int use_sparse;
extern int transient_method;
static double *dc_point;
extern struct components_t *g_components;
extern struct instruction_t *g_instructions;
extern struct option_t *g_options;

void instruction_tran_tr(struct instruction_t *instr, int max_nodes, int sources, 
    double *MNA_G, double *MNA_C, double *RHS, double *L, 
    double *U, double *m, int *P,double *temp)
{
  struct instruction_t *ptr;
  double t;
  int i,j;

  double *matrix = (double*) z_malloc(sizeof(double)*(max_nodes+sources)*(max_nodes+sources));
  double *matrix2 = (double*) z_malloc(sizeof(double)*(max_nodes+sources)*(max_nodes+sources));
  double *solution_0 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *solution_1 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *RHS_0 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *RHS_1 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *swap;


  
  if ( P )
    for (i=0; i<max_nodes+sources; i++ )
      for ( j=0; j<max_nodes+sources; j++ )
        temp[i*(max_nodes+sources)+j] = MNA_G[ P[i]*(max_nodes+sources)+j];
  
  swap = temp;
  temp = MNA_G;
  MNA_G = swap;

  memcpy(RHS_0, RHS, sizeof(double)*(max_nodes+sources));
  memcpy(solution_0, dc_point, sizeof(double)*(max_nodes+sources));

  matrix_multiply_scalar(MNA_C, MNA_C, 2/instr->tran.time_step, max_nodes+sources);
  matrix_add_matrix(matrix, MNA_G, MNA_C, (max_nodes+sources));
  matrix_sub_matrix(matrix2, MNA_G, MNA_C, (max_nodes+sources));
  
  
  if ( P )
    for( i = 0 ; i < max_nodes + sources; i++ )
      P[i] = i;

  if ( iter_type == NoIter ) {
    memset(L, 0, sizeof(double)*(max_nodes+sources));
    memset(U, 0, sizeof(double)*(max_nodes+sources));

    if ( spd_flag==0) 
      LU_decomposition(matrix, L, U, P, max_nodes + sources );
    else {
      Choleski_LDU_Decomposition(matrix, L, max_nodes + sources );
      calculate_transpose(L, U, max_nodes + sources );
    }
  } else {
    for (i=0; i<max_nodes+sources; i++ ) 
      m[i] = matrix[i*(max_nodes+sources)+i];
  }

  for ( t=instr->tran.time_step; t<=instr->tran.time_finish; t+= instr->tran.time_step ) {
    calculate_RHS(g_components, max_nodes, sources, RHS_1, t, 0);
    multiply_matrix_vector(matrix2,solution_0, temp, max_nodes+sources);
    add_vectors(RHS_0, RHS_1, RHS_0, max_nodes+sources);
    sub_vectors(RHS_0, temp, RHS_0, max_nodes+sources);
    printf("---\n");
    print_array(RHS_0, max_nodes+sources);
    if ( iter_type == NoIter ) {
      forward_substitution(L, RHS_0, temp ,P, max_nodes + sources);
      backward_substitution(U, temp, solution_1, max_nodes+sources);
    } else {
      memset(solution_1, 0, sizeof(double)*(max_nodes+sources));
      if ( iter_type==CG)
        conjugate(matrix, solution_1, RHS_0, m, itol, max_nodes+sources);
      else
        biconjugate(matrix, solution_1, RHS_0, m, itol, max_nodes+sources);
    }
    
    printf("RHS: \n");

    for (i=0; i< max_nodes + sources ;i++ )	
      printf("%7g\n", RHS_1[i]);

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
  
  free(matrix);
  free(matrix2);
  free(solution_0);
  free(solution_1);
  free(RHS_0);
  free(RHS_1);
}

void instruction_tran_be(struct instruction_t *instr, int max_nodes, int sources, 
    double *MNA_G, double *MNA_C, double *RHS, double *L, 
    double *U, double *m, int *P,double *temp)
{
  struct instruction_t *ptr;
  double t;
  int i,j;

  double *matrix = (double*) z_malloc(sizeof(double)*(max_nodes+sources)*(max_nodes+sources));
  double *solution_0 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *solution_1 = (double*) z_malloc(sizeof(double)*(max_nodes+sources));
  double *swap;
  
  if ( P )
    for (i=0; i<max_nodes+sources; i++ )
      for ( j=0; j<max_nodes+sources; j++ )
        temp[i*(max_nodes+sources)+j] = MNA_G[ P[i]*(max_nodes+sources)+j];
  
  swap = temp;
  temp = MNA_G;
  MNA_G = swap;

  memcpy(solution_0, dc_point, sizeof(double)*(max_nodes+sources));

  matrix_multiply_scalar(MNA_C, MNA_C, 1/instr->tran.time_step, max_nodes+sources);
  matrix_add_matrix(matrix, MNA_G, MNA_C, (max_nodes+sources));
    
  if ( P )
    for( i = 0 ; i < max_nodes + sources; i++ )
      P[i] = i;

  if ( iter_type == NoIter ) {
    if ( spd_flag==0) 
      LU_decomposition(matrix, L, U, P, max_nodes + sources );
    else {
      Choleski_LDU_Decomposition(matrix, L, max_nodes + sources );
      calculate_transpose(L, U, max_nodes + sources );
    }
  } else {
    for (i=0; i<max_nodes+sources; i++ ) 
      m[i] = matrix[i*(max_nodes+sources)+i];
  }

  for ( t=instr->tran.time_step; t<=instr->tran.time_finish; t+= instr->tran.time_step ) {
    calculate_RHS(g_components, max_nodes, sources, RHS, t, 0);
    multiply_matrix_vector(MNA_C,solution_0, temp, max_nodes+sources);
    add_vectors(RHS, temp, RHS, max_nodes+sources);

    if ( iter_type == NoIter ) {
      forward_substitution(L, RHS, temp ,P, max_nodes + sources);
      backward_substitution(U, temp, solution_1, max_nodes+sources);
    } else {
      memset(solution_1, 0, sizeof(double)*(max_nodes+sources));
      if ( iter_type==CG)
        conjugate(matrix, solution_1, RHS, m, itol, max_nodes+sources);
      else
        biconjugate(matrix, solution_1, RHS, m, itol, max_nodes+sources);
    }

    /*printf("RHS: \n");

    for (i=0; i< max_nodes + sources ;i++ )	
      printf("%7g\n", RHS_0[i]);*/

    /*printf("Result:\n");
    print_array(solution_1,max_nodes+sources);*/

    for ( ptr = g_instructions; ptr!=NULL; ptr=ptr->next) {

      if ( ptr->type == Plot ) {
        for (j=0; j<ptr->plot.num; j++ )
          fprintf( ptr->plot.output[j], "%7G %g\n", 
              t, solution_1[ptr->plot.list[j]]);

      }
    }

    swap = solution_0;
    solution_0 = solution_1;
    solution_1 = swap;
  }

  free(matrix);
  free(solution_0);
  free(solution_1);

}

void instruction_tran(struct instruction_t *instr, int max_nodes, int sources, 
    int renamed_nodes[], double *MNA_G, double *MNA_C, double *RHS, double *L, 
    double *U, double *m, int *P,double *temp)
{
  if ( transient_method == TR ) {
    instruction_tran_tr(instr, max_nodes, sources, MNA_G, MNA_C, RHS, L, U, m, P, temp);
  } else if ( transient_method == BE ) {
    instruction_tran_be(instr, max_nodes, sources, MNA_G, MNA_C, RHS, L, U, m, P, temp);
  } else {
    printf("Unknown transient method :%d\n", transient_method);
  }
}


void circuit_mna(struct components_t *circuit, double **MNA_G, double **MNA_C, int *max_nodes, int *sources,
    int element_types[], int **renamed_nodes)
{
  int transient_analysis = 0;
  struct components_t *s;
  struct instruction_t *w;
  int max_v_id, x, y;
  int  elements;
  double *matrix, *matrix2;

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
        if ( s->data.t1.is_ground==0 )
          (*sources)++;
        break;
      case R:
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

  for ( w = g_instructions; w; w=w->next ){
    if ( w->type == Tran ) {
      transient_analysis = 1;
      break;
    }
  }


  *MNA_G= (double*) z_calloc((*max_nodes+*sources)*(*max_nodes+*sources), sizeof(double));


  if ( transient_analysis )
    *MNA_C= (double*) z_calloc((*max_nodes+*sources)*(*max_nodes+*sources), sizeof(double));

  matrix  = *MNA_G;
  matrix2 = *MNA_C;

  // meta ftiaxnoume to panw aristera elements x elements pou einai ta pa8htika stoixeia
  for (s=circuit; s!=NULL ;s=s->next) {
    if ( s->data.type == R ) {
      if ( s->data.t1.plus < *max_nodes)
        matrix[ s->data.t1.plus + s->data.t1.plus*(*max_nodes+*sources)] += 1/s->data.t1.val;
      if ( s->data.t1.minus < *max_nodes )
        matrix[ s->data.t1.minus + s->data.t1.minus*(*max_nodes+*sources)] += 1/s->data.t1.val;
      if ( s->data.t1.minus < *max_nodes && s->data.t1.plus < *max_nodes ) {
        matrix[ s->data.t1.minus + s->data.t1.plus*(*max_nodes+*sources)] -= 1/s->data.t1.val;
        matrix[ s->data.t1.plus + s->data.t1.minus*(*max_nodes+*sources)] -= 1/s->data.t1.val;
      }
    } else if ( transient_analysis && s->data.type == C )  {
      if ( s->data.t1.plus < *max_nodes)
        matrix2[ s->data.t1.plus + s->data.t1.plus*(*max_nodes+*sources)] += s->data.t1.val;
      if ( s->data.t1.minus < *max_nodes )
        matrix2[ s->data.t1.minus + s->data.t1.minus*(*max_nodes+*sources)] += s->data.t1.val;
      if ( s->data.t1.minus < *max_nodes && s->data.t1.plus < *max_nodes ) {
        matrix2[ s->data.t1.minus + s->data.t1.plus*(*max_nodes+*sources)] -= s->data.t1.val;
        matrix2[ s->data.t1.plus + s->data.t1.minus*(*max_nodes+*sources)] -= s->data.t1.val;
      }
    } else if ( transient_analysis && s->data.type == L ) {
      matrix2[*max_nodes+*sources-inductors +s->data.t1.id 
        + (*max_nodes+*sources)*(*max_nodes+*sources-inductors+s->data.t1.id)] = - s->data.t1.val;
    }

    if ( ( s->data.type == V && s->data.t1.is_ground == 0 )  || s->data.type == L ) {
      if ( s->data.t1.plus < *max_nodes ) {
        matrix[ *max_nodes + s->data.t1.id + (*max_nodes+*sources) * s->data.t1.plus ] = 1;
        matrix[ ( *max_nodes + s->data.t1.id ) * (*max_nodes+*sources) + s->data.t1.plus ] = 1;
      }
      if ( s->data.t1.minus < *max_nodes ) {
        matrix[ *max_nodes + s->data.t1.id + (*max_nodes+*sources) * s->data.t1.minus] = -1;
        matrix[ ( *max_nodes + s->data.t1.id ) * (*max_nodes+*sources) + s->data.t1.minus] = -1;
      }
    }
  }
  
  FILE* mna = fopen("mna_analysis", "w");

  printf("MNA matrix: see file \"mna_analysis\"\n" );

  for (y= 0; y < *max_nodes+*sources; y++ ) {
    /*if ( y == *max_nodes )
      fprintf(mna, "\n");*/
    for ( x=0; x < *max_nodes+*sources; x++ ) {
      /*if ( x == *max_nodes )
        fprintf(mna, "  ");*/
      fprintf(mna, "%10g", matrix[x + y*(*max_nodes+*sources)] );

    }
    fprintf(mna, "\n");
  }
  fclose(mna);

  if ( transient_analysis) { 
    mna = fopen("mna_analysis_transient", "w");

    printf("MNA matrix: see file \"mna_analysis_transient\"\n" );

    for (y= 0; y < *max_nodes+*sources; y++ ) {
      /*if ( y == *max_nodes )
        fprintf(mna, "\n");*/
      for ( x=0; x < *max_nodes+*sources; x++ ) {
        /*if ( x == *max_nodes )
          fprintf(mna, "  ");*/
        fprintf(mna, "%10g", matrix2[x + y*(*max_nodes+*sources)] );

      }
      fprintf(mna, "\n");
    }
    fclose(mna);

  }
}


void solve(double *L, double *U, double *temp, double *result,
    double *RHS,
    int *P, int max_nodes, int sources, double t, int dc_only)
{
  calculate_RHS(g_components,max_nodes,sources,RHS, t, dc_only);
  forward_substitution(L, RHS, temp ,P, max_nodes + sources);
  backward_substitution(U, temp, result, max_nodes+sources);
}


int instruction_dc(struct instruction_t *instr, int max_nodes, int sources, int renamed_nodes[], 
    double *MNA, double *RHS, double *L, double *U, double *m, int *P,
    double *temp, double *result)
{
  double begin;
  double end;
  double step;
  double dummy;
  int i;
  int original_source_val;
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
          printf("Solving for Voltage Source value %g\n",dummy);
          if ( iter_type == NoIter )
            solve(L,U,temp,result,RHS,P,max_nodes,sources, -1, 1);
          else {
            calculate_RHS(g_components,max_nodes,sources,RHS, -1, 1);
            memset(result, 0, sizeof(double) * ( max_nodes+sources));
            if ( iter_type == CG ) {
              conjugate(MNA, result , RHS, m, itol, max_nodes+sources);
            } else if ( iter_type == BiCG ) {
              biconjugate(MNA, result , RHS, m, itol, max_nodes+sources);
            }
          }
          printf("Result:\n");
          print_array(result,max_nodes+sources);

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
          printf("Solve for Current Source value %g\n",dummy);


          if ( iter_type == NoIter ) {
            solve(L,U,temp,result,RHS,P,max_nodes,sources,-1, 1);
          } else {
            calculate_RHS(g_components,max_nodes,sources,RHS, -1, 1);
            memset(result, 0, sizeof(double) * ( max_nodes+sources));
            if ( iter_type == CG ) {
              conjugate(MNA, result , RHS, m, itol, max_nodes+sources);
            } else if ( iter_type == BiCG ) {
              biconjugate(MNA, result , RHS, m, itol, max_nodes+sources);
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

  return 0;
}

int execute_instructions(double *MNA_G, double *MNA_C,  int max_nodes, int sources, int *renamed_nodes, int stoixeia[])
{
  double *RHS = NULL;
  double *L=NULL,*U=NULL,*result=NULL, *m=NULL, *temp=NULL;
  int *P=NULL;
  int i;
  struct instruction_t *instr;
  double *test = (double*) calloc(max_nodes+sources,sizeof(double));

  // Gia tous sparse pinakes

  instr = g_instructions;

  RHS = (double*) z_calloc(max_nodes+sources, sizeof(double));
  result = (double*) z_calloc((max_nodes+sources), sizeof(double));
  dc_point = (double*) z_calloc((max_nodes+sources), sizeof(double));
  temp = (double* ) z_calloc( (max_nodes+sources) * (max_nodes+sources), sizeof(double));

  if ( iter_type == NoIter ) {
    L = (double*) z_calloc((max_nodes+sources)*(max_nodes+sources), sizeof(double));
    U = (double*) z_calloc((max_nodes+sources)*(max_nodes+sources), sizeof(double));
    P = (int*) z_calloc((max_nodes+sources), sizeof(int));
    for( i = 0 ; i < max_nodes + sources; i++ )
      P[i] = i;
    if ( spd_flag==0) {
      LU_decomposition(MNA_G, L, U, P, max_nodes + sources );
    } else {
      if ( stoixeia[typeV] ) {
        printf("[-] Den prepei na uparxoun phges Tashs\n");
        return 0;
      }

      if ( stoixeia[typeL] ) {
        printf("[-] Den prepei na uparxoun phnia\n");
        return 0;
      }

      if ( stoixeia[typeD] ) {
        printf("[-] Den prepei na uparxoun diodoi\n");
        return 0;
      }

      if ( stoixeia[typeQ] ) {
        printf("[-] Den prepei na uparxoun transistor BJT \n");
        return 0;
      }

      if ( stoixeia[typeM] ) {
        printf("[-] Den prepei na uparxoun transistor CMOS\n");
        return 0;
      }


      if ( !Choleski_LDU_Decomposition(MNA_G, L, max_nodes + sources ) ) {
        printf("[-] Negative value on array L\n");
        exit(0);
      }

      calculate_transpose(L, U, max_nodes + sources );
    }

    solve(L,U,temp,dc_point,RHS,P,max_nodes,sources, -1, 1);
    printf("Circuit Solution\n");
    print_array(dc_point, max_nodes+sources);

  } else {
    m = (double*) z_malloc(sizeof(double) * (max_nodes+sources));

    for (i=0; i<max_nodes+sources; i++ ) 
      m[i] = MNA_G[i*(max_nodes+sources)+i];

    calculate_RHS(g_components,max_nodes,sources,RHS, -1, 1);
    biconjugate(MNA_G, dc_point, RHS, m, itol, max_nodes+sources);
    printf("Circuit Solution\n");
    print_array(dc_point, max_nodes+sources);
  }
calculate_RHS(g_components,max_nodes,sources,RHS, -1, 1);

  FILE *dc = fopen("dc_point", "w");
  if ( P ) {
    for (i=0; i<max_nodes+sources; i++)
      fprintf(dc,"%10g\n", dc_point[P[i]]);
    fprintf(dc,"...........RHS...\n");
    for (i=0; i<max_nodes+sources; i++)
      fprintf(dc,"%10g\n", RHS[P[i]]);
  } else {
    for (i=0; i<max_nodes+sources; i++)
      fprintf(dc,"%10g\n", dc_point[i]);
    fprintf(dc,"...........RHS...\n");
    for (i=0; i<max_nodes+sources; i++)
      fprintf(dc,"%10g\n", RHS[i]);
  }
    

  multiply_matrix_vector(MNA_G,dc_point,test,max_nodes+sources);

  fprintf(dc,"...........\n");

  for(i = 0 ; i < max_nodes+sources ; i++)
    fprintf(dc,"%10g\n",test[i]);

  fclose(dc);
  free(test);


  while ( instr ) {
    switch ( instr->type ) {
      case Dc:
        instruction_dc(instr,max_nodes, sources, renamed_nodes, 
            MNA_G, RHS, L, U, m, P, temp, result); 
        break;

      case Tran:
        instruction_tran(instr, max_nodes, sources, renamed_nodes,
            MNA_G, MNA_C, RHS, L, U, m, P, temp);
        break;
      case Plot:
      break;
    }
    instr = instr->next;
  }

  if (dc_point )
    free(dc_point);

  if ( m ) free(m);
  if ( L ) free(L);
  if ( U ) free(U);
  if ( P ) free(P);
  if ( temp ) free(temp);
  if ( result) free(result);
  if ( MNA_G) free(MNA_G);
  if ( MNA_C) free(MNA_C);
  if ( RHS ) free(RHS);
  if ( renamed_nodes ) free(renamed_nodes);

  return 0;
}
