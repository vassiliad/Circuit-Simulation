%{
#include <stdio.h>
#include <stdlib.h>
#include "types.h"
#include "syntactic.tab.h"

#define YYERROR_VERBOSE 1
extern int yylineno;
void yyerror(const char *str);
int yylex();
void new_t1(struct T1_t* ret, int id, int plus, int minus, double val);
void new_t2(struct T2_t* ret, int id, int plus, int minus, char *model_name, int area_used, double area);
void new_t3(struct T3_t* ret, int id, int d, int g, int s, int b, double l, double w, char* model_name);
void new_t4(struct T4_t* ret, int id, int c, int b, int e, char *model_name, int area_used, double area);

double number_to_double(const number_t *number);


struct components_t *g_components;
struct instruction_t *g_instructions;
struct option_t *g_options;

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


%token V I R C L D Q M NUMBER NEW_LINE STR
%token V2
%token PLOT DC OPTION TRAN
%token ASSIGN

%token EXP SIN PWL PULSE
%token LPAREN RPAREN

%type <str> STR
%type <description> description component
%type <tail1> tail1
%type <tail2> tail2
%type <tail3> tail3
%type <tail4> tail4
%type <number> NUMBER;
%type <id>V I R C L D Q M V2 
%type <entries> entries;
%type <instruction> instruction
%type <num_list> v_sources_list;
%type <option> option;
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
}
;

entries: entries option {
	$2->prev = g_options;
	if ( g_options )
		g_options -> next = $2;
	g_options = $2;
	$$ = $1;
}
| entries instruction {
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
| instruction {
	$$.instructions = $1;
	g_instructions = $1;
}
| component
{
	g_components = (struct components_t*) calloc(1,sizeof(struct components_t));
	g_components->data = $1;

	$$.components = g_components;
}
| option {
	g_options = $1;
}

option: OPTION STR {
	$$ = (struct option_t * ) calloc(1,sizeof(struct option_t));

	if ( strcasecmp($2,"spd") == 0 ) {
		$$->type= SPD;
	} else if ( strcasecmp($2, "iter") == 0 ) {
    $$->type = ITER;
    $$->iter_type = CG;
  } else {
		yyerror("[-] Unknown option\n");
	}
  free($2);
}
| OPTION STR STR
{
	$$ = (struct option_t * ) calloc(1,sizeof(struct option_t));
  if (strcasecmp($2, "iter")==0) {
    $$->type = ITER;
    if ( strcasecmp($3, "spd") == 0 )
      $$->iter_type = BiCG;
    else {
      yyerror("[-] Invalid iter type\n");
    }
  }
  free($2);
  free($3);
}
| OPTION STR ASSIGN NUMBER
{
	$$ = (struct option_t * ) calloc(1,sizeof(struct option_t));
  if ( strcasecmp($2, "itol") == 0 ) {
    $$->type = ITOL;
    $$->itol = ( $4.type == Integer ? $4.integer : $4.dbl );
  } else {
    yyerror("[-] Unknown Option\n");
		free($$);
		return 1;
  }

  free($2);
}
| OPTION STR ASSIGN STR
{
	$$ = ( struct option_t * ) calloc(1, sizeof(struct option_t));
	
	if ( strcasecmp($2, "method") == 0 ) {
		if ( strcasecmp($4, "tr") == 0 ) {
			$$->type = TR;
		} else if ( strcasecmp($4, "be" ) == 0 ) {
			$$->type = BE;
		} else {
			yyerror("[-] Expected \"TR\" or \"BE\"\n");
			free($$);
			return 1;
		}
		free($2);
		free($4);
	} else {
		yyerror("[-] Expected \"METHOD\"\n");
		free($$);
		return 1;
	}
}
;


instruction: DC V NUMBER NUMBER NUMBER {
	$$ = (struct instruction_t*) calloc(1,sizeof(struct instruction_t));
	$$->type = Dc;
	$$->dc.sourceType = Voltage;
	$$->dc.source = $2;

	$$->dc.begin = $3.type == Integer ? $3.integer : $3.dbl;
	$$->dc.end   = $4.type == Integer ? $4.integer : $4.dbl;
	$$->dc.inc   = $5.type == Integer ? $5.integer : $5.dbl;
}
| DC I NUMBER NUMBER NUMBER {
	$$ = (struct instruction_t*) calloc(1,sizeof(struct instruction_t));
	$$->type = Dc;
	$$->dc.sourceType = Current;
	$$->dc.source = $2;

	$$->dc.begin = $3.type == Integer ? $3.integer : $3.dbl;
	$$->dc.end   = $4.type == Integer ? $4.integer : $4.dbl;
	$$->dc.inc   = $5.type == Integer ? $5.integer : $5.dbl;
	
	$$->prev = NULL;

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
			printf("[-] Could not open %s for output (plot)\n", filename);
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
	$$.vals = (int*) realloc($$.vals, sizeof(int) * ($$.num+1) );
	$$.vals[ $$.num ] = $2;
	$$.num++;
}
| V2 {
	$$.vals = (int*) calloc(1,sizeof(int));
	$$.vals[0] = $1;
	$$.num = 1;
}
;

component: description NEW_LINE { $$ = $1; }
;

description:V tail1 transient_spec 
{
	$$.type = V;
	new_t1(&($$.t1),$1, $2.plus, $2.minus, $2.val);

	if ( $2.val == 0 )
		$$.t1.is_ground = 1;
	
	$$.t1.transient = $3;
}
| I tail1 transient_spec {
	$$.type = I;
	new_t1(&($$.t1),$1, $2.plus, $2.minus, $2.val);

	$$.t1.transient = $3;
// void new_t1(struct T1_t* ret, int id, int plus, int minus, double val)
}
| R tail1 {
	$$.type = R;
	new_t1(&($$.t1),$1, $2.plus, $2.minus, $2.val);
// void new_t1(struct T1_t* ret, int id, int plus, int minus, double val)
}
| C tail1 {
	$$.type = C;
	new_t1(&($$.t1),$1, $2.plus, $2.minus, $2.val);
// void new_t1(struct T1_t* ret, int id, int plus, int minus, double val)
}
| L tail1 {
	$$.type = L;
	new_t1(&($$.t1),$1, $2.plus, $2.minus, $2.val);
// void new_t1(struct T1_t* ret, int id, int plus, int minus, double val)
}
| D tail2 {
	$$.type = D;
	new_t2(&($$.t2), $1, $2.plus, $2.minus,  $2.model_name, 0, 0  );
//void new_t2(struct T2_t* ret, int id, int plus, int minus, char *model_name, int area_used, double area)
}
| D tail2 NUMBER {
	if ( $3.type != Double )
		yyerror("Expected Double as last argument\n");
	$$.type = D;
	new_t2(&($$.t2), $1, $2.plus, $2.minus, $2.model_name,  1, $3.dbl );
//void new_t2(struct T2_t* ret, int id, int plus, int minus, char *model_name, int area_used, double area)
}
| M tail3 {
	$$.type = M;
	new_t3(&($$.t3), $1, $2.d, $2.g, $2.b, $2.s, $2.l, $2.w, $2.model_name);
//void new_t3(struct T3_t* ret, int id, int d, int g, int s, int b, double l, double w, char* model_name)
}
| Q tail4 {
	$$.type = Q;
	new_t4(&($$.t4), $1, $2.c, $2.b, $2.e, $2.model_name, 0, 0 );
	//void new_t4(struct T4_t* ret, int id, int c, int b, int e, char *model_name, int area_used, double area)
}
| Q tail4 NUMBER{
	if ( $3.type != Double )
			yyerror("Expected Double as last argument\n");
	$$.type = Q;
	new_t4(&($$.t4), $1, $2.c, $2.b, $2.e, $2.model_name, 1, $3.dbl );
	//void new_t4(struct T4_t* ret, int id, int c, int b, int e, char *model_name, int area_used, double area)
}
;

tail4: NUMBER NUMBER NUMBER STR {
	if ( $1.type != Integer )
		yyerror("C should be integer\n");
	if ( $2.type != Integer )
		yyerror("B should be integer\n");
	if ( $3.type != Integer )
		yyerror("E should be integer\n");
	$$.c = $1.integer;
	$$.b = $2.integer;
	$$.e = $3.integer;
	$$.model_name  = $4;
}
;

tail3: NUMBER NUMBER NUMBER NUMBER   NUMBER NUMBER STR {
	//INT INT INT INT DOUBLE DOUBLE STR

	if ( $1.type != Integer )
		yyerror("D should be integer\n");
	if ( $2.type != Integer)
		yyerror("G should be integer\n");
	if ( $3.type != Integer)
		yyerror("S should be integer\n");
	if ( $4.type != Integer )
		yyerror("B should be integer\n");
	
	if ( $5.type != Double )
		yyerror("L should be double\n");
	if ($6.type != Double )
		yyerror("W should be double\n");

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
	if ($1.type != Integer )
		yyerror("PLUS should be integer\n");
	if ($2.type != Integer )
		yyerror("MINUS should be integer\n");
	
	$$.plus = $1.integer;
	$$.minus = $2.integer;
	$$.model_name = $3;
}
;

tail1: NUMBER NUMBER NUMBER {
	if ($1.type != Integer )
		yyerror("PLUS should be integer\n");
	if ($2.type != Integer )
		yyerror("MINUS should be integer\n");
	
	$$.plus = $1.integer;
	$$.minus = $2.integer;
	
	if ( $3.type == Integer )
		$$.val = $3.integer;
	else
		$$.val = $3.dbl;
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
	ret->is_ground = 0; // auto einai gia tis phges taseis se ola ta alla
											// den exei shmasia h timh tou
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
	fprintf(stderr,"Error: %s@%d\n", str,yylineno);
}

