%option noyywrap
%option yylineno
%{
#include <stdio.h>
#include <ctype.h>
#include "types.h"
#include "syntactic.tab.h"

#define zprint(...)

%}

%option nounput
%option noinput

%%

\*.*[\n]? { }

"," return COMMA;
"(" return LPAREN;
")" return RPAREN;

"=" {
  return ASSIGN;
}

".OPTION" {
	return OPTION;
}

".TRAN" return TRAN;

".PLOT" {
	return PLOT;
}

".DC" {
	return DC;
}

".op" {
}

".end" {
}

"PWL" return PWL;
"SIN" return SIN;
"PULSE" return PULSE;
"EXP" return EXP;

[-+]?[[:digit:]]+ {
	yylval.number.type = Integer;
	yylval.number.integer = strtoul(yytext, 0, 10);
  zprint("Number_i: %d\n", yylval.number.integer);
	return NUMBER;
}

[-+]?[[:digit:]]+(\.[[:digit:]]+(e[-+]?[[:digit:]]+)?|(\.[[:digit:]]+)?(e[-+]?[[:digit:]]+)) {
	yylval.number.type = Double;
	yylval.number.dbl = strtod(yytext,NULL);
  zprint("Number_d: %g\n", yylval.number.dbl);
	return NUMBER;
}



V[ \t]*\([ \t]*[[:digit:]]+[ \t]*\) {
  char * str;
  char *ptr;

  str = strchr(yytext, '(')+1;
  ptr = strchr(yytext, ')');
  str = strndup(str, ptr-str);
  
  zprint("V2 %s (%s)\n", yytext, str);
  yylval.str = str;
	return V2;
}

\n+ {
	return NEW_LINE;
}

[\ \r\t]+ { }

[_]?[[:alpha:]]+([[:alnum:]]?[_]?)+ {
	yylval.str = strdup(yytext);
  zprint("str: %s\n", yytext);
	return STR;
}

. { printf("Unexpected character %s@%d\n", yytext,yylineno); exit(0);}
%%

