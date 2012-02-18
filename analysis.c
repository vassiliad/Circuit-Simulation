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
extern double *dc_point;
extern struct components_t *g_components;
extern struct instruction_t *g_instructions;
extern struct option_t *g_options;



void circuit_mna(struct components_t *circuit, double **MNA_G, double **MNA_C, int *max_nodes, int *sources,
    int element_types[], int **renamed_nodes)
{
  int transient_analysis = 0;
  struct components_t *s, *p;
  struct instruction_t *w;
  int max_v_id;
  int  elements;
  int x,y,c;
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
        if ( s->data.t1.val > 0 )
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


  *MNA_G= (double*) calloc((*max_nodes+*sources)*(*max_nodes+*sources), sizeof(double));


  if ( transient_analysis )
    *MNA_C= (double*) calloc((*max_nodes+*sources)*(*max_nodes+*sources), sizeof(double));

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

    if ( s->data.type == V  || s->data.type == L ) {
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
    int *P, int max_nodes, int sources, double t)
{
  calculate_RHS(g_components,max_nodes,sources,RHS, t);
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
            solve(L,U,temp,result,RHS,P,max_nodes,sources, -1);
          else {
            calculate_RHS(g_components,max_nodes,sources,RHS, -1);
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
            solve(L,U,temp,result,RHS,P,max_nodes,sources,-1);
          } else {
            calculate_RHS(g_components,max_nodes,sources,RHS, -1);
            memset(result, 0, sizeof(double) * ( max_nodes+sources));
            if ( iter_type == CG ) {
              conjugate(MNA, result , RHS, m, itol, max_nodes+sources);
            } else if ( iter_type == BiCG ) {
              biconjugate(MNA, result , RHS, m, itol, max_nodes+sources);
            }
          }
          printf("Result:\n");
          print_array(result,max_nodes+sources);

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

  // Gia tous sparse pinakes

  instr = g_instructions;

  RHS = (double*) calloc(max_nodes+sources, sizeof(double));
  result = (double*) calloc((max_nodes+sources), sizeof(double));
  dc_point = (double*) calloc((max_nodes+sources), sizeof(double));

  if ( iter_type == NoIter ) {
    L = (double*) calloc((max_nodes+sources)*(max_nodes+sources), sizeof(double));
    U = (double*) calloc((max_nodes+sources)*(max_nodes+sources), sizeof(double));
    P = (int*) calloc((max_nodes+sources), sizeof(int));
    temp = (double* ) calloc( (max_nodes+sources) * (max_nodes+sources), sizeof(double));
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

    solve(L,U,temp,dc_point,RHS,P,max_nodes,sources, -1);
    printf("Circuit Solution\n");
    print_array(dc_point, max_nodes+sources);

  } else {
    m = (double*) malloc(sizeof(double) * (max_nodes+sources));

    for (i=0; i<max_nodes+sources; i++ ) 
      m[i] = MNA_G[i*(max_nodes+sources)+i];

    calculate_RHS(g_components,max_nodes,sources,RHS, -1);
    biconjugate(MNA_G, dc_point, RHS, m, itol, max_nodes+sources);
    printf("Circuit Solution\n");
    print_array(dc_point, max_nodes+sources);
  }
  


  while ( instr ) {
    if(instr->type == Dc) {
      instruction_dc(instr,max_nodes, sources, renamed_nodes, 
	  MNA_G, RHS, L, U, m, P, temp, result); 
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
  if ( RHS ) free(RHS);
  if ( renamed_nodes ) free(renamed_nodes);

  return 0;
}
