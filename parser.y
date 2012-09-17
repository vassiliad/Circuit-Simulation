%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "options.h"
#include "components.h"
#include "transient.h"
#include "hash_table.h"
#include "plot.h"

#define IDS_CHUNK 1000

#define YYERROR_VERBOSE 1
extern int yylineno;
int yylex();
int yyerror(const char *str);

%}

%union
{
  int integer;
  double dbl;
  char *string;
  transient_t *transient;
  pair_t pair;
  pwl_t pwl;
};
%error-verbose
%token NEW_LINE COMMA INTEGER DOUBLE STRING PLOT_V
%token DC OPTIONS TRAN PLOT
%token EXP SIN PWL PULSE
%token LPAREN RPAREN ASSIGN


%type <integer> INTEGER node_id
%type <dbl> DOUBLE number
%type <string> STRING PLOT_V
%type <transient> transient_spec
%type <pwl> pairs 
%type <pair> pair


%%

spice: entries
{
};

number: DOUBLE
{
  $$ = $1;
}
| INTEGER {
  $$ = (double) $1;
}
;

node_id: STRING
{
  $$ = hash_get($1);
}
| INTEGER
{
  if ( $1 == 0 ) {
    $$ = 0;
  } else {
    char temp[50];
    sprintf(temp, "%u", $1);
    $$ = hash_get(strdup(temp));
  }
}
;

entries: entries component NEW_LINE
| entries options NEW_LINE
| entries instruction NEW_LINE
| entries NEW_LINE
| component NEW_LINE
| options NEW_LINE
| NEW_LINE
;

component: 
STRING node_id node_id number transient_spec
{
  switch(tolower($1[0])) {
    case 'v' :
      new_v( $2, $3, $4, $5);
    break;

    case 'i':
      new_i( $2, $3, $4, $5);
    break;

    case 'r':
      if ( $5 )
        return yyerror("Cannot have transient_spec in a resistor");
      new_r( $2, $3, $4);
    break;

    case 'c':
      if ( $5 )
        return yyerror("Cannot have transient_spec in a capasitor");
      new_c( $2, $3, $4);
    break;

    case 'l':
      if ( $5 )
        return yyerror("Cannot have transient_spec in an inductor");
      new_l( $2, $3, $4);
    break;

    default:
      free($1);
      return yyerror("Unknown component");
  }
  free($1);
}

;

options: OPTIONS option_list
{
}

plot_list: plot_list PLOT_V
{
  char *c = strdup($2);
  plot_node(hash_get($2), c);
  free(c);
}
| PLOT_V
{
  char *c = strdup($1);
  plot_node(hash_get($1), c);
  free(c);
}
;

instruction: TRAN number number
{
  do_transient = 1;
  tran_step = $2;
  tran_finish = $3;
}
| PLOT plot_list
;

pairs: pairs pair {
  $$.pairs = (pair_t*) realloc($$.pairs, sizeof(pair_t)*($$.size+1));
  $$.pairs[$$.size] = $2;
  $$.size++;
}
| pair
{
  $$.pairs = (pair_t*) calloc(1, sizeof(pair_t));
  *($$.pairs)=$1;
  $$.size = 1;
}

pair: LPAREN number number RPAREN
{
  $$.t = $2;
  $$.i = $3;
}

transient_spec: SIN LPAREN number number number number number number RPAREN
{
  $$ = (transient_t*) calloc(1, sizeof(transient_t));
  $$->type = Sin;
  $$->tsin.i1 = $3;
  $$->tsin.ia = $4;
  $$->tsin.fr = $5;
  $$->tsin.td = $6;
  $$->tsin.df = $7;
  $$->tsin.ph = $8;
}
| EXP LPAREN number number number number number number RPAREN
{
  $$ = ( transient_t*) calloc(1, sizeof(transient_t));

  $$->type = Exp;

  $$->texp.i1 = $3;
  $$->texp.i2 = $4;
  $$->texp.td1 = $5;
  $$->texp.tc1 = $6;
  $$->texp.td2 = $7;
  $$->texp.tc2 = $8;
}
| PULSE LPAREN number number number number number number number RPAREN
{
  $$ = ( transient_t* ) calloc(1, sizeof(transient_t));
  $$->type = Pulse;

  $$->tpulse.i1 = $3;
  $$->tpulse.i2 = $4;
  $$->tpulse.td = $5;
  $$->tpulse.tr = $6;
  $$->tpulse.tf = $7;
  $$->tpulse.pw = $8;
  $$->tpulse.per = $9;
}
| PWL pairs
{
   $$ = ( transient_t* ) calloc(1, sizeof(transient_t));
   $$->type = Pwl;
   
   $$->tpwl = $2;
}
|
{
  $$ = NULL;
}
;

option_list: option_list COMMA option
| option
;

option: STRING
{
  if ( strcasecmp($1, "spd") == 0 ) {
    method_noniter = CholDecomp;
  } else if (strcasecmp($1, "sparse")==0) {
    sparse_use = 1;
  } else if (strcasecmp($1, "iter") == 0) {
    method_choice = Iterative;
    method_iter = CG;
  } else{
    free($1);
    return yyerror("Uknown option");
  }
  free($1);
}
| STRING STRING
{
  if (strcasecmp($1, "iter")==0) {
    method_choice = Iterative;
    method_iter   = CG;
    if ( strcasecmp($2, "spd") == 0 )
      method_iter = BiCG;
    else {
      yyerror("Invalid iter type");
      free($1);
      free($2);
      return 1;
    }
  }
  free($1);
  free($2);
}
| STRING ASSIGN STRING
{
  if ( strcasecmp($1, "method") == 0 ) {
		if ( strcasecmp($3, "tr") == 0 ) {
			method_tran = Tr;
		} else if ( strcasecmp($3, "be" ) == 0 ) {
			method_tran = Be;
		} else {
			yyerror("Expected \"TR\" or \"BE\"");
			free($1);
      free($3);
			return 1;
		}
	} else {
		yyerror("Expected \"METHOD\"");
    free($1);
    free($3);
		return 1;
	}
  free($1);
  free($3);
}
| STRING ASSIGN number
{
  if ( strcasecmp($1, "itol") == 0 ) {
    itol = $3;
  } else {
    yyerror("Unknown Option");
    free($1);
		return 1;
  }

  free($1);
}
;
%%

int yyerror(const char *str)
{
  printf("Error[%d]: %s\n", yylineno, str);
  return 1;
}
