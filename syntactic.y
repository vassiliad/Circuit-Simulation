%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "types.h"
#include "syntactic.tab.h"

#define IDS_CHUNK 1000

#define YYERROR_VERBOSE 1
extern int yylineno;
void yyerror(const char *str);
int yylex();
void new_t1(struct T1_t* ret, int id, int plus, int minus, double val);
void new_t2(struct T2_t* ret, int id, int plus, int minus, char *model_name, int area_used, double area);
void new_t3(struct T3_t* ret, int id, int d, int g, int s, int b, double l, double w, char* model_name);
void new_t4(struct T4_t* ret, int id, int c, int b, int e, char *model_name, int area_used, double area);

double number_to_double(const number_t *number);
int compare_pair(const void *p1, const void *p2);

struct components_t *g_components;
struct instruction_t *g_instructions;
struct option_t *g_options;

int str_to_id(char *str);
void free_ids();

int num_ids = 0;
int max_num_ids = 0;
char **ids;


%}

%union 
{
	int id;					// V C I R L Q D M
	struct NUMBER_T number;
	char *str;

	struct T1_t t1; // tail1 V C I R L
	struct T2_t t2;	// tail2 D
	struct T3_t t3; // tail3 M
	struct T4_t t4; // tail4 Q

	struct description_t description;

	struct components_t *components;
	
	struct tail1_t {	// V C I R L
		int plus, minus;
		double val;
	}tail1;

	struct tail2_t {	// D
		int plus, minus;	
		char* model_name;
	} tail2;

	struct tail3_t {	// Q
		int s,g,b,d;
		double l,w;
		char *model_name;
	}tail3;

	struct tail4_t { // M
		int c,b,e;
		char* model_name;
	} tail4;
	
	struct {
		int *vals;
		int num;
	} num_list;

	entries_t entries;
	struct instruction_t *instruction;
	struct option_t *option;
	transient_spec_t *transient;
	pwl_t pairs;
	pair_t pair;
};


%token NUMBER NEW_LINE STR COMMA
%token V2
%token PLOT DC OPTION TRAN
%token ASSIGN

%token EXP SIN PWL PULSE
%token LPAREN RPAREN

%type <str> STR V2
%type <description> description component
%type <tail1> tail1
%type <tail2> tail2
%type <tail3> tail3
%type <tail4> tail4
%type <number> NUMBER;
%type <entries> entries;
%type <instruction> instruction
%type <num_list> v_sources_list
%type <option> option options option_entry
%type <transient> transient_spec
%type <pairs> pairs
%type <pair> pair
%start input_file

%error-verbose
%%
input_file: entries {
	struct components_t *s = g_components;
	struct instruction_t *p =g_instructions;
	struct option_t *c = g_options;
  struct components_t *ground;

  ground = ( struct components_t * ) calloc(1, sizeof(struct components_t));
  // add ground
  ground->data.type = V;
	new_t1(&(ground->data.t1),-1, 0, 0,0);
  ground->data.t1.is_ground = 1;

  if ( s ) 
    s ->next = ground;
  else {
    g_components = s = ground;
  }

	/* rewind the list, components is the list's tail */
	while ( s && s->prev )
		s = s->prev;
	while (p && p->prev)
		p= p->prev;
	while ( c && c->prev )
		c = c->prev;

	g_options = c;
	g_instructions = p;
	g_components = s;

  free_ids();
}
;

entries: entries option NEW_LINE {
	$2->prev = g_options;
	if ( g_options )
		g_options -> next = $2;
	g_options = $2;
	$$ = $1;
}
| entries instruction NEW_LINE {
	$2->prev = g_instructions;
	if ( g_instructions )
		g_instructions->next = $2;
	g_instructions = $2;
	$$ = $1;
	$$.instructions = g_instructions;
}
| entries NEW_LINE{
	$$ = $1;
}
| entries component {
	struct components_t *comp = (struct components_t*) calloc(1,sizeof(struct components_t));
	comp->data = $2;

	comp->prev = g_components;
	if ( g_components )
		g_components->next = comp;
	
	g_components = comp;
	$$ = $1;
	$$.components = g_components;
}
| instruction  NEW_LINE{
	$$.instructions = $1;
	g_instructions = $1;
}
| component
{
	g_components = (struct components_t*) calloc(1,sizeof(struct components_t));
	g_components->data = $1;

	$$.components = g_components;
}
| option NEW_LINE{
	g_options = $1;
}

option: OPTION options
{
  $$ = $2;
}

options: options COMMA option_entry
{
  struct option_t *o;

  for ( o= $1; o->next; o = o->next );
  
  o->next = $3;
  $$ = $1;
}
| option_entry
{
  $$ = $1;
}
option_entry: STR  {
	$$ = (struct option_t * ) calloc(1,sizeof(struct option_t));

	if ( strcasecmp($1,"spd") == 0 ) {
		$$->type= SPD;
	}  else if (strcasecmp($1, "sparse") == 0 ) {
    $$->type = SPARSE;
  } else if ( strcasecmp($1, "iter") == 0 ) {
    $$->type = ITER;
    $$->iter_type = CG;
  } else {
		yyerror("Unknown option");
    return 1;
	}
  free($1);
}
| STR STR
{
	$$ = (struct option_t * ) calloc(1,sizeof(struct option_t));
  if (strcasecmp($1, "iter")==0) {
    $$->type = ITER;
    if ( strcasecmp($2, "spd") == 0 )
      $$->iter_type = BiCG;
    else {
      yyerror("Invalid iter type");
      return 1;
    }
  }
  free($1);
  free($2);
}
| STR ASSIGN NUMBER
{
	$$ = (struct option_t * ) calloc(1,sizeof(struct option_t));
  if ( strcasecmp($1, "itol") == 0 ) {
    $$->type = ITOL;
    $$->itol = ( $3.type == Integer ? $3.integer : $3.dbl );
  } else {
    yyerror("Unknown Option");
    free($1);
		free($$);
		return 1;
  }

  free($1);
}
| STR ASSIGN STR
{
	$$ = ( struct option_t * ) calloc(1, sizeof(struct option_t));
	
	if ( strcasecmp($1, "method") == 0 ) {
		if ( strcasecmp($3, "tr") == 0 ) {
			$$->type = TR;
		} else if ( strcasecmp($3, "be" ) == 0 ) {
			$$->type = BE;
		} else {
			yyerror("Expected \"TR\" or \"BE\"");
			free($$);
			return 1;
		}
		free($1);
		free($3);
	} else {
		yyerror("Expected \"METHOD\"");
		free($$);
    free($1);
    free($3);
		return 1;
	}
}
;


instruction: DC STR NUMBER NUMBER NUMBER {
  // V I
	$$ = (struct instruction_t*) calloc(1,sizeof(struct instruction_t));
	$$->type = Dc;
  if ( tolower(($2)[0]) == 'v' )
  	$$->dc.sourceType = Voltage;
  else if ( tolower(($2)[0]) =='i' )
    $$->dc.sourceType = Current;
  else {
    yyerror("DC is only for V/I");
    return 1;
  }

	$$->dc.source = str_to_id(strdup(($2)+1));

	$$->dc.begin = $3.type == Integer ? $3.integer : $3.dbl;
	$$->dc.end   = $4.type == Integer ? $4.integer : $4.dbl;
	$$->dc.inc   = $5.type == Integer ? $5.integer : $5.dbl;
}
| PLOT v_sources_list {
	int i = 0;
	char filename[256];
	$$ = (struct instruction_t*) calloc(1,sizeof(struct instruction_t));
	$$->type = Plot;
	$$->plot.list = $2.vals;
	$$->plot.num = $2.num;
	$$->plot.output = (FILE**) calloc(1,sizeof(FILE*) * $$->plot.num);

	for (i=0; i<$$->plot.num; i++ ) {
		sprintf(filename, "plot_v_%d", $2.vals[i]);
		$$->plot.output[i] = fopen(filename, "w");
		if ( $$->plot.output[i] == NULL ) {
			printf("Could not open %s for output (plot)\n", filename);
			exit(0);
		}
	}
}
| TRAN NUMBER NUMBER
{
	$$ = ( struct instruction_t *) calloc(1, sizeof(struct instruction_t));
	$$->type = Tran;
	$$->tran.time_step = number_to_double(&$2);
	$$->tran.time_finish = number_to_double(&$3);
}
;

v_sources_list: v_sources_list V2 {
  char *endptr;

	$$.vals = (int*) realloc($$.vals, sizeof(int) * ($$.num+1) );
	$$.vals[$$.num] = strtol($2, &endptr, 10);

  if ( *endptr!=0 && !isspace(*endptr) ) {
  	$$.vals[$$.num] = str_to_id($2);
  }
  $$.num++;
}
| V2 {
  char *endptr;
	$$.vals = (int*) calloc(1,sizeof(int));
  
  $$.vals[0] = strtol($1, &endptr, 10);

  if ( *endptr!=0 && !isspace(*endptr) ) {
  	$$.vals[0] = str_to_id($1);
  }
	$$.num = 1;
}
;

component: description NEW_LINE { $$ = $1; }
;

description:STR tail1 transient_spec 
{
  // V I
	
  switch(tolower($1[0])) {
    case 'v':
      $$.type = V;
    break;

    case 'i':
      $$.type = I;
    break;

    case 'r':
      if ( $3 != 0 ) {
        yyerror("R cannot have transient_spec\n");
        return 1;
      }
      $$.type = R;
    break;
    
    case 'l':
      if ( $3 != 0 ) {
        yyerror("L cannot have transient_spec\n");
        return 1;
      }
      $$.type = L;
    break;
    case 'c':
      if ( $3 != 0 ) {
        yyerror("C cannot have transient_spec\n");
        return 1;
      }
      $$.type = C;
    break;
    
    default: {
      char err[500];
      sprintf(err, "Invalid component %c", (char)toupper(($1)[0]));
      yyerror(err);
      return 1;
    }

  }
	new_t1(&($$.t1),str_to_id(strdup($1+1)), $2.plus, $2.minus, $2.val);
  free($1);
	$$.t1.transient = $3;
}
| STR tail2 {
  // d
	$$.type = D;
	new_t2(&($$.t2), str_to_id(strdup($1+1)), $2.plus, $2.minus,  $2.model_name, 0, 0  );
  free($1);
}
| STR tail2 NUMBER {
  // d
	if ( $3.type != Double ) {
		yyerror("Expected Double as last argument");
    return 1;
  }
	$$.type = D;
	new_t2(&($$.t2), str_to_id(strdup($1+1)), $2.plus, $2.minus, $2.model_name,  1, $3.dbl );
  free($1);
}
| STR tail3 {
  // M
	$$.type = M;
	new_t3(&($$.t3), str_to_id(strdup($1+1)), $2.d, $2.g, $2.b, $2.s, $2.l, $2.w, $2.model_name);
  free($1);
//void new_t3(struct T3_t* ret, int id, int d, int g, int s, int b, double l, double w, char* model_name)
}
| STR tail4 {
  // Q
	$$.type = Q;
	new_t4(&($$.t4), str_to_id(strdup($1+1)), $2.c, $2.b, $2.e, $2.model_name, 0, 0 );
  free($1);
	//void new_t4(struct T4_t* ret, int id, int c, int b, int e, char *model_name, int area_used, double area)
}
| STR tail4 NUMBER{
  // Q
	if ( $3.type != Double ) {
			yyerror("Expected Double as last argument"); 
      return 1;
  }
	$$.type = Q;
	new_t4(&($$.t4), str_to_id(strdup($1+1)), $2.c, $2.b, $2.e, $2.model_name, 1, $3.dbl );
  free($1);
	//void new_t4(struct T4_t* ret, int id, int c, int b, int e, char *model_name, int area_used, double area)
}
;

tail4: NUMBER NUMBER NUMBER STR {
	if ( $1.type != Integer ) {
		yyerror("C should be integer");
    return 1;
  }
	if ( $2.type != Integer ) {
		yyerror("B should be integer");
    return 1;
  }
	if ( $3.type != Integer ) {
		yyerror("E should be integer");
    return 1;
  }
	$$.c = $1.integer;
	$$.b = $2.integer;
	$$.e = $3.integer;
	$$.model_name  = $4;
}
;

tail3: NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER STR {
	//INT INT INT INT DOUBLE DOUBLE STR

	if ( $1.type != Integer ) {
		yyerror("D should be integer");
    return 1;
  }

	if ( $2.type != Integer) {
		yyerror("G should be integer");
    return 1;
  }

	if ( $3.type != Integer) {
		yyerror("S should be integer");
    return 1;
  }

	if ( $4.type != Integer ) {
		yyerror("B should be integer");
    return 1;
  }
	
	if ( $5.type != Double ) {
		yyerror("L should be double");
    return 1;
  }

	if ($6.type != Double ) {
		yyerror("W should be double");
    return 1;
  }

	$$.d = $1.integer;
	$$.g = $2.integer;
	$$.s = $3.integer;
	$$.b = $4.integer;
	$$.l = $5.dbl;
	$$.w = $6.dbl;
	$$.model_name = $7;
}
;

tail2: NUMBER NUMBER STR {
	if ($1.type != Integer ) {
		yyerror("PLUS should be integer");
    return 1;
  }

	if ($2.type != Integer ) {
		yyerror("MINUS should be integer");
    return 1;
  }
	
	$$.plus = $1.integer;
	$$.minus = $2.integer;
	$$.model_name = $3;
}
;


tail1: NUMBER NUMBER NUMBER {
  if ( $1.type != Integer) {
    yyerror("Plus should be an integer");
    return 1;
  }
  if ( $2.type != Integer ) {
    yyerror("Minus should be an integer");
    return 1;
  }
  $$.plus = $1.integer;

  $$.minus = $2.integer;
  $$.val = number_to_double(&$3);
} 
| STR NUMBER NUMBER
{ 
  if ( $2.type != Integer ) {
    yyerror("Minus should be an integer");
    return 1;
  }

  $$.plus = str_to_id($1);
  $$.minus = $2.integer;
  $$.val = number_to_double(&$3);
}
| NUMBER STR NUMBER
{
  if ( $1.type != Integer) {
    yyerror("Plus should be an integer");
    return 1;
  }
  $$.plus = $1.integer;
  $$.minus = str_to_id($2);
  $$.val = number_to_double(&$3);
}
| STR STR NUMBER
{
  $$.plus = str_to_id($1);
  $$.minus = str_to_id($2);
  $$.val = number_to_double(&$3);
}
;

transient_spec: EXP LPAREN NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER RPAREN
{
	$$ = ( transient_spec_t* ) calloc(1, sizeof(transient_spec_t));
	$$->type = Exp;
	$$->exp.i1 = number_to_double(&$3);
	$$->exp.i2 = number_to_double(&$4);
	$$->exp.td1 = number_to_double(&$5);
	$$->exp.td2 = number_to_double(&$6);
	$$->exp.tc1 = number_to_double(&$7);
	$$->exp.tc2 = number_to_double(&$8);
}
| SIN LPAREN NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER RPAREN
{
  $$ = ( transient_spec_t* ) calloc(1, sizeof(transient_spec_t));
	$$->type = Sin;
	$$->_sin.i1 = number_to_double(&$3);
	$$->_sin.ia = number_to_double(&$4);
	$$->_sin.fr = number_to_double(&$5);
	$$->_sin.td = number_to_double(&$6);
	$$->_sin.df = number_to_double(&$7);
	$$->_sin.ph = number_to_double(&$8);
}
| PULSE LPAREN NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER RPAREN
{
	$$ = ( transient_spec_t* ) calloc(1, sizeof(transient_spec_t));
	$$->type = Pulse;
	$$->pulse.i1 = number_to_double(&$3);
	$$->pulse.i2 = number_to_double(&$4);
	$$->pulse.td = number_to_double(&$5);
	$$->pulse.tr = number_to_double(&$6);
	$$->pulse.tf = number_to_double(&$7);
	$$->pulse.pw = number_to_double(&$8);
	$$->pulse.per = number_to_double(&$9);
}
| PWL pairs
{
	$$ = ( transient_spec_t* ) calloc(1, sizeof(transient_spec_t));
	$$->type = Pwl;
	$$->pwl.pairs = $2.pairs;
	$$->pwl.size  = $2.size;
  
  if ( $2.size < 2 ) {
    yyerror("PWL Pairs must be at least 2");
    return 1;
  }
  qsort($$->pwl.pairs,$2.size, sizeof(pair_t), compare_pair);
}
|
{
	$$ = NULL;
};

pairs: pairs pair
{
	$$ = $1;
	$$.size ++;
	$$.pairs = ( pair_t* ) realloc($$.pairs, sizeof(pair_t)*($$.size));
	$$.pairs[ $$.size-1 ] = $2;
}
| pair
{
	$$.size = 1 ;
	$$.pairs = ( pair_t * ) calloc(1, sizeof(pair_t));
	$$.pairs[0] = $1;
}

pair: LPAREN NUMBER NUMBER RPAREN
{
	$$.i = number_to_double(&$2);
	$$.t = number_to_double(&$3);
}

%%

double number_to_double(const number_t *number)
{
	if ( number->type == Integer )
		return number->integer;
	else
		return number->dbl;
}

/*****Helper func, to init the struct*****/
void new_t1(struct T1_t* ret, int id, int plus, int minus, double val)
{
	ret->transient = NULL;
	ret->id = id;
	ret->plus = plus;
	ret->minus = minus;
	ret->val = val;
  ret->is_ground = 0;
}

/*****Helper func, to init the struct*****/
void new_t2(struct T2_t* ret, int id, int plus, int minus, char *model_name, int area_used, double area)
{
	ret->id = id;
	ret->plus = plus;
	ret->minus = minus;
	ret->model_name = model_name;
	ret->area_used = area_used;
	ret->area = area;
}

/*****Helper func, to init the struct*****/
void new_t3(struct T3_t* ret, int id, int d, int g, int s, int b, double l, double w, char* model_name)
{
	ret->id = id;
	ret->d = d;
	ret->g = g;
	ret->s = s;
	ret->b = b;
	ret->l = l;
	ret->w = w;
	ret->model_name = model_name;
}

/*****Helper func, to init the struct*****/
void new_t4(struct T4_t* ret, int id, int c, int b, int e, char *model_name, int area_used, double area)
{
	ret->id=id;
	ret->c=c;
	ret->b=b;
	ret->e=e;
	ret->model_name=model_name;
	ret->area_used=area_used;
	ret->area=area;
}

void yyerror(const char *str)
{
	fprintf(stderr,"[-] Error[%d]: %s\n",yylineno, str);
}

int compare_pair(const void *p1, const void *p2)
{
  if ( ((pair_t*)p1)->t > ((pair_t*)p2)->t )
    return 1;
  if ( ((pair_t*)p1)->t < ((pair_t*)p2)->t )
    return -1;

  return 0;
}


int str_to_id(char *str)
{
  int i =0;
  int j;
  // strip whitespace from str

  for (i=0; str[i]!=0; i++ ) {
    if ( isspace(str[i]) ) {
      for (j=i; str[j] != 0; j++)
        str[j] = str[j+1];
      i--;
    }
  }

  //printf("str: %s\n", str);

  for ( i=0 ;i<num_ids; i++ )
    if ( strcasecmp(ids[i], str) == 0 )
      return i;

  if ( num_ids == max_num_ids ) {
    max_num_ids+= IDS_CHUNK;
  }
  ids = (char**) realloc(ids, sizeof(char*) * (max_num_ids));
  ids[num_ids] = str;

  return num_ids++;
}

void free_ids()
{
  int i;
  for (i=0; i<num_ids; i++ )
    free(ids[i]);
  free(ids);

  num_ids = 0;
}

