#include "utility.h"
#include "types.h"
extern int spd_flag;
extern int iter_type;
extern double itol;
extern int use_sparse;
extern struct components_t *g_components;
extern struct instruction_t *g_instructions;
extern struct option_t *g_options;




double linear_interpolate(double x, double x0, double y0, double x1, double y1)
{
  return y0 + ((x-x0)*y1 - (x-x0)*y0)/(x1-x0);
}

double calculate_ac(const transient_spec_t *transient, double t)
{
  int i;
  if ( transient == NULL ) 
    return 0;

  switch ( transient->type ) {
    case Exp:
      if ( t < transient->exp.td1 ) {
        return transient->exp.i1;
      } else if ( t < transient->exp.td2 ) {
        return transient->exp.i1 + ( transient->exp.i2 - transient->exp.i1 )*
            ( 1 - pow(epsilon, -(t-transient->exp.td1)/transient->exp.tc1));
      } else {
        return transient->exp.i1 + ( transient->exp.i2 - transient->exp.i1 )*
            ( pow(epsilon, -(t-transient->exp.td2)/transient->exp.tc2) 
            - pow(epsilon, -(t-transient->exp.td1)/transient->exp.tc1));
      }
    break;

    case Sin:
      if ( t < transient->_sin.td ) {
        return transient->_sin.i1 + transient->_sin.ia * sin( 2*pi * transient->_sin.ph / 360);
      } else {
        return transient->_sin.i1 + transient->_sin.ia 
          * sin(2*pi * transient->_sin.fr * (t - transient->_sin.td ) 
          + 2*pi * transient->_sin.ph / 360) 
          * pow(epsilon, -(t-transient->_sin.td)* transient->_sin.df);
      }

    break;
    case Pulse:
      t-= transient->pulse.per * (int)( t/transient->pulse.per);

      if ( t < transient->pulse.td ) {
        return transient->pulse.i1;
      } else if ( t  < transient->pulse.td + transient->pulse.tr ) {
        return  linear_interpolate(t,
                                  transient->pulse.td,
                                  transient->pulse.i1,
                                  transient->pulse.td+transient->pulse.tr,
                                  transient->pulse.i2);
      } else if ( t  < transient->pulse.td + transient->pulse.tr + 
                       transient->pulse.pw) {
        return transient->pulse.i2;
      } else if ( t  < transient->pulse.td + transient->pulse.tr +
                       transient->pulse.pw + transient->pulse.tf ) {
        return linear_interpolate(t,
                                  transient->pulse.td+transient->pulse.tr+transient->pulse.pw,
                                  transient->pulse.i2,
                                  transient->pulse.td+transient->pulse.tr+transient->pulse.pw
                                      + transient->pulse.tf,
                                  transient->pulse.i1);
      } else {
        return transient->pulse.i1;
      }
       
    break;

    case Pwl:
      for ( i=0; i<transient->pwl.size; i++ ) {
        if ( t > transient->pwl.pairs[i].t ) {
          if( i == transient->pwl.size-1 ) {
            return transient->pwl.pairs[i].i;            
          }

          return linear_interpolate(t,
                      transient->pwl.pairs[i].t,
                      transient->pwl.pairs[i].i,
                      transient->pwl.pairs[i+1].t,
                      transient->pwl.pairs[i+1].i);
        }
      }

      return transient->pwl.pairs->i;
    break;

    default:
      printf("Invalid calculate_ac type\n");
      exit(0);
  }
  
}


void options_cleanup(struct option_t *g_options)
{
  struct option_t *o = g_options, *n;

  while ( o ) {
    n = o->next;
    free(o);
    o = n;
  }
}

void instructions_cleanup(struct instruction_t *instr)
{
  struct instruction_t  *next;
  int i =0;

  while ( instr != NULL ) {
    next = instr->next;

    if ( instr->type == Plot ) {
      free( instr->plot.list );
      for (i=0; i< instr->plot.num; i++ )
        fclose(instr->plot.output[i]);
      free( instr->plot.output );
    }

    free(instr);

    instr = next;
  }
}

void circuit_print(struct components_t *circuit)
{
  struct components_t *s;
  s = g_components;

  for (;s!=NULL; s=s->next) {
    switch ( s->data.type ) {
      case V:
        printf("V"); 
        printf("%d +:%d -:%d v:%lf\tis ground:%s\n", s->data.t1.id, s->data.t1.plus, s->data.t1.minus, s->data.t1.val,
            (s->data.t1.val==0 ? "yes" : "no"));
        break;
      case I:
        printf("I");
        printf("%d +:%d -:%d i:%lf\n", s->data.t1.id, s->data.t1.plus, s->data.t1.minus, s->data.t1.val);
        break;
      case R:
        printf("R");
        printf("%d +:%d -:%d r:%lf\n", s->data.t1.id, s->data.t1.plus, s->data.t1.minus, s->data.t1.val);
        break;
      case C:
        printf("C");
        printf("%d +:%d -:%d c:%lf\n", s->data.t1.id, s->data.t1.plus, s->data.t1.minus, s->data.t1.val);
        break;
      case L:
        printf("L");
        printf("%d +:%d -:%d l:%lf\n", s->data.t1.id, s->data.t1.plus, s->data.t1.minus, s->data.t1.val);
        break;
      case D:
        printf("D%d +:%d -:%d model:%s", s->data.t2.id, s->data.t2.plus, s->data.t2.minus, s->data.t2.model_name);
        if ( s->data.t2.area_used )
          printf(" area:%lf\n", s->data.t2.area );
        else
          printf("\n");
        break;
      case Q:
        printf("Q%d c:%d b:%d e:%d model:%s", s->data.t4.id, s->data.t4.c, s->data.t4.b, s->data.t4.e, s->data.t4.model_name);
        if ( s->data.t4.area_used )
          printf(" area:%lf\n", s->data.t4.area );
        else
          printf("\n");
        break;
      case M:
        printf("M%d d:%d g:%d s:%d b:%d l:%lf w:%lf model:%s\n", s->data.t3.id, s->data.t3.d, s->data.t3.g, 
            s->data.t3.s, s->data.t3.b, s->data.t3.l, s->data.t3.w, s->data.t3.model_name);
        break;
      default: printf("unknown: %d !!!\n", s->data.type);
    }
  }
}


int calculate_RHS(struct components_t *circuit,int max_nodes,int sources, double *RHS, double t){
  int y;
  struct components_t *s;
  double temp = 0;

  for ( y=0; y< max_nodes + sources; y++ )
    RHS[y] = 0;

  for (s=circuit; s!=NULL; s=s->next) {
    if ( s->data.type == I ) {
      if ( t >= 0 )
        temp = calculate_ac(s->data.t1.transient, t);

      if ( s->data.t1.plus != max_nodes )
        RHS[ s->data.t1.plus ] += s->data.t1.val + temp;
      else if ( s->data.t1.minus != max_nodes )
        RHS[ s->data.t1.minus] -= s->data.t1.val - temp;
      else {
        printf("Vraxukuklwma sthn phgh reumatos panw sth geiwsh\n");
        exit(0);
      }
      

    }
    else if ( s->data.type == V && s->data.t1.is_ground == 0 ) {
      if ( t >= 0 )
        temp = calculate_ac(s->data.t1.transient, t);
      RHS[ max_nodes + s->data.t1.id ] = s->data.t1.val + temp;
    }
  }

    return 0;
}

void print_help(char *path)
{
  char *file;
  file = strrchr(path, '/');

  if ( file == NULL )
    file = path;
  else
    file ++;

  printf("Usage: %s SOURCE_FILE\n", file);
}



void instructions_print(struct instruction_t *instr)
{
  int i;
  if ( instr == NULL ) {
    printf("[#] No instructions\n");
    return;
  }
  while ( instr ) {
    switch(instr->type) {
      case Dc:
        printf("DC source: %d[%s] begin: %g end: %g step: %g\n",
            instr->dc.source, (instr->dc.sourceType == Voltage? "Voltage" : "Current" ), instr->dc.begin, instr->dc.end, instr->dc.inc);
        break;

      case Plot:
        printf("PLOT sources: " );
        for ( i=0; i<instr->plot.num-1; i++ )
          printf("%d, ", instr->plot.list[i]);
        printf("%d\n", instr->plot.list[i]);
        break;

      default:
        printf("Unsupported instruction (%d)\n", instr->type);
    }
    instr = instr->next;
  }
}





void circuit_cleanup(struct components_t *circuit)
{
  struct components_t *e, *p;

  e = circuit;

  while ( e ) {
    switch ( e->data.type ) {
      case D:
        free( e->data.t2.model_name);
        break;
      case M:
        free(e->data.t3.model_name);
        break;
      case Q:
        free(e->data.t4.model_name);
        break;
    }
    p = e->next;
    free(e);
    e = p;
  }

}






int circuit_rename_nodes(struct components_t *circuit, struct instruction_t *instr, int element_types[], int **renamed_nodes, int *max_nodes)
{
  int *indices = NULL;
  int *voltage_sources = NULL, num_voltage_sources=0;
  int *current_sources = NULL, num_current_sources=0;
  int num_indices = 0, i,j;
  struct components_t *s;

  element_types[typeV] = element_types[typeI] = element_types[typeR] = element_types[typeC]
    = element_types[typeL] = element_types[typeD] = element_types[typeQ] = element_types[typeM]
    = 0;


  struct components_t *ground = NULL;
  struct instruction_t *p;


  for (s=circuit; s!=NULL; s=s->next) {
    switch( s->data.type ) {
      case V:
        if ( s->data.t1.is_ground )
          ground = s;
        else {
          voltage_sources = (int*) realloc(voltage_sources, sizeof(int)*(num_voltage_sources+1) );
          voltage_sources[ num_voltage_sources++ ] = s->data.t1.id;
          s->data.t1.original_id = s->data.t1.id;
          s->data.t1.id = element_types[typeV]++;
        }
        break;
      case I:
        current_sources =(int*)realloc(current_sources, sizeof(int)* ( num_current_sources+1));
        current_sources[ num_current_sources++ ] = s->data.t1.id;
        s->data.t1.original_id = s->data.t1.id;
        s->data.t1.id = element_types[typeI]++;
        break;
      case R:
        s->data.t1.id = element_types[typeR]++;
        break;
      case L:
        s->data.t1.id = element_types[typeL]++;
        break;
      case C:
        s->data.t1.id = element_types[typeC]++;
        break;
      case D:
        s->data.t2.id = element_types[typeD]++;
        break;
      case M:
        s->data.t3.id = element_types[typeM]++;
        break;
      case Q:
        s->data.t3.id = element_types[typeQ]++;
        break;
      default:
        printf("Unknown element type in rename (%d)\n", s->data.type);
        continue;
    }

    switch ( s->data.type ) {
      case V: case I: case R: case L: case C:
        for (i=0; i<num_indices; i++)
          if ( indices[i] == s->data.t1.plus ) {
            s->data.t1.plus = i;
            break;
          }

        if ( i == num_indices ) { // paei na pei oti den to brhka mesa sto indices
          indices = (int*) realloc(indices, sizeof(int) * (++num_indices));
          indices[num_indices-1] = s->data.t1.plus;
          s->data.t1.plus = num_indices-1;
        }

        for (i=0; i<num_indices; i++)
          if ( indices[i] == s->data.t1.minus) {
            s->data.t1.minus = i;
            break;
          }

        if ( i == num_indices ) { // paei na pei oti den to brhka mesa sto indices
          indices = (int*) realloc(indices, sizeof(int) * (++num_indices));
          indices[num_indices-1] = s->data.t1.minus;
          s->data.t1.minus = num_indices-1;
        }
        break;

      default:
        printf("[#] Den exei ginei akoma gia ta upoloipa sto indices mesa sthn rename\n");
        return 0;
        break;
    }
  }

  if ( ground == NULL ) 
  {
    printf("To do8en kuklwma den exei geiwsh\n");
    exit(0);
  }

  if ( ground != NULL ) {
    ground->data.t1.id = element_types[typeV];
    int ground_node = ground->data.t1.plus;

    int temp = indices[ground_node];
    indices[ground_node] = indices[num_indices-1];
    indices[num_indices-1] = temp;

    for (s=circuit; s!=NULL; s=s->next) {
      switch( s->data.type ) {
        case V:
        case I:
        case R:
        case C:
        case L:
          if ( s->data.t1.plus == ground_node )
            s->data.t1.plus = num_indices-1;
          else if ( s->data.t1.plus == num_indices-1 )
            s->data.t1.plus = ground_node;
          if ( s->data.t1.minus == ground_node)
            s->data.t1.minus = num_indices-1;
          else if ( s->data.t1.minus == num_indices-1 )
            s->data.t1.minus = ground_node;
          break;
        case M:
          if ( s->data.t3.d == ground_node )
            s->data.t3.d = num_indices-1;
          else if ( s->data.t3.d == num_indices-1)
            s->data.t3.d = ground_node;

          if ( s->data.t3.g == ground_node )
            s->data.t3.g = num_indices-1;
          else if ( s->data.t3.g == num_indices-1)
            s->data.t3.g = ground_node;

          if ( s->data.t3.s == ground_node )
            s->data.t3.s = num_indices-1;
          else if ( s->data.t3.s == num_indices-1)
            s->data.t3.s = ground_node;

          if ( s->data.t3.b == ground_node )
            s->data.t3.b = num_indices-1;
          else if ( s->data.t3.b == num_indices-1)
            s->data.t3.b = ground_node;

          break;

        case Q:
          if ( s->data.t4.c == ground_node )
            s->data.t4.c = num_indices-1;
          else if ( s->data.t4.c == num_indices-1)
            s->data.t4.c = ground_node;

          if ( s->data.t4.b == ground_node )
            s->data.t4.b = num_indices-1;
          else if ( s->data.t4.b == num_indices-1)
            s->data.t4.b = ground_node;

          if ( s->data.t4.e == ground_node )
            s->data.t4.e = num_indices-1;
          else if ( s->data.t4.e == num_indices-1)
            s->data.t4.e = ground_node;
          break;

        case D:
          if ( s->data.t2.plus == ground_node )
            s->data.t2.plus = num_indices-1;
          else if ( s->data.t2.plus == num_indices-1 )
            s->data.t2.plus = ground_node;
          if ( s->data.t2.minus == ground_node)
            s->data.t2.minus = num_indices-1;
          else if ( s->data.t2.minus == num_indices-1 )
            s->data.t2.minus = ground_node;

          break;
        default:
          printf("[-] Error sto rename den uposthrizetai o tupos akoma! (%d)\n", s->data.type);
          break;
      }
    }
    *max_nodes = num_indices-1;
    *renamed_nodes = indices;
  } else {
    *max_nodes = num_indices;
    *renamed_nodes = indices;
  }

  printf("Nodes that were renamed:\n");
  for (i=0; i<num_indices; i++ )
    printf("\t%d -> %d\n", indices[i], i);

  for (p=instr; p!=NULL; p=p->next) {
    switch(p->type) {
      case Dc:
        if ( p->dc.sourceType == Voltage ) {
          for (i=0; i<num_voltage_sources; i++ ) {
            if ( voltage_sources[i] == p->dc.source ) {
              p->dc.source = i;
              break;
            }
          }

          if ( i == num_voltage_sources ) { // auto paei na pei oti den bre8hke h phgh
            printf("[-] Could not find voltage source %d (for instruction DC)\n", p->dc.source);
            exit(0);
          }
        } else if ( p->dc.sourceType == Current ) {
          for (i=0; i<num_current_sources; i++ ) {
            if ( current_sources[i] == p->dc.source ) {
              p->dc.source = i;
              break;
            }
          }

          if ( i == num_current_sources ) { // auto paei na pei oti den bre8hke h phgh
            printf("[-] Could not find current source %d (for instruction DC)\n", p->dc.source);
            exit(0);
          }
        }

        break;

      case Plot:
        for ( i=0; i< p->plot.num; i++ ) {
          for (j=0; j<num_indices; j++ ) {
            if ( indices[j] == p->plot.list[i] ) {
              if ( ground && j == num_indices-1 ) {
                printf("[#] Will not report Ground's voltage\n");
                fclose(p->plot.output[i]);

                for (j=i; j<p->plot.num-1; j++ ) {
                  p->plot.list[j] = p->plot.list[j+1];
                  p->plot.output[j] = p->plot.output[j+1];
                }
                p->plot.num--;
                p->plot.list = (int*) realloc(p->plot.list, sizeof(int) * ( p->plot.num));
                p->plot.output = (FILE**) realloc(p->plot.output, sizeof(FILE*) * ( p->plot.num));
                i--;
              } else {
                printf("plot ( was %d turned to %d )\n", p->plot.list[i], j);
                p->plot.list[i] = j;
              }
              break;
            }
          }

          if ( j == num_indices ) {

            printf("[-] Plot cannot report %d's voltage because it does not exist\n", p->plot.list[i]);
            exit(0);
          }
        }
        break;

      default:
        printf("[#] Rename Nodes --> Rename VoltageSources in instructions ( UNKNOWN INSTRUCTION %d )\n", p->type);
        break;
    }
  }

  if ( voltage_sources)
    free(voltage_sources);
  if (current_sources)
    free(current_sources);

  if ( ground )
    return num_indices-1;
  else
    return num_indices;
}
