%option noyywrap
%option yylineno
%{
#include <stdio.h>
#include <ctype.h>
#include "types.h"
#include "syntactic.tab.h"

%}
%%

\*.*[\n]? { }

"(" return LPAREN;
")" return RPAREN;

"=" {
  return ASSIGN;
}

".OPTION" {
	return OPTION;
}

".PLOT" {
	return PLOT;
}

".DC" {
	return DC;
}

"PWL" return PWL;
"SIN" return SIN;
"PULSE" return PULSE;
"EXP" return EXP;

[[:digit:]]+ {
	yylval.number.type = Integer;
	yylval.number.integer = strtoul(yytext, 0, 10);
	return NUMBER;
}

[[:digit:]]+(\.[[:digit:]]+(e[-+]?[[:digit:]]+)?|(\.[[:digit:]]+)?(e[-+]?[[:digit:]]+)) {
	yylval.number.type = Double;
	yylval.number.dbl = strtod(yytext,NULL);
	return NUMBER;
}



V[ \t]*[[:digit:]]+ {
	char *s = yytext;
	while ( isdigit(*s) == 0 )
		s++;
	yylval.id = strtoul(s,0,10);
	return V;
}

I[ \t]*[[:digit:]]+  {
	char *s = yytext;
	while ( isdigit(*s) == 0 )
		s++;
	yylval.id = strtoul(s,0,10);
	return I;
}

V[ \t]*\([ \t]*[[:digit:]]+[ \t]*\) {
	char *s = yytext;
	while ( isdigit(*s) == 0 )
		s++;
	yylval.id = strtoul(s,0,10);
	return V2;
}

R[ \t]*[[:digit:]]+  {
	char *s = yytext;
	while ( isdigit(*s) == 0 )
		s++;
	yylval.id = strtoul(s,0,10);
	return R;
}

C[ \t]*[[:digit:]]+  {
	char *s = yytext;
	while ( isdigit(*s) == 0 )
		s++;
	yylval.id = strtoul(s,0,10);
	return C;
}

L[ \t]*[[:digit:]]+  {
	char *s = yytext;
	while ( isdigit(*s) == 0 )
		s++;
	yylval.id = strtoul(s,0,10);
	return L;
}

D[ \t]*[[:digit:]]+  {
	char *s = yytext;
	while ( isdigit(*s) == 0 )
		s++;
	yylval.id = strtoul(s,0,10);
	return D;
}

Q[ \t]*[[:digit:]]+  {
	char *s = yytext;
	while ( isdigit(*s) == 0 )
		s++;
	yylval.id = strtoul(s,0,10);
	return Q;
}

M[ \t]*[[:digit:]]+  {
	char *s = yytext;
	while ( isdigit(*s) == 0 )
		s++;
	yylval.id = strtoul(s,0,10);
	return M;
}


\n+ {
	return NEW_LINE;
}

" " { }
[\ \r\t]+ { }

[_]?[[:alpha:]]+([[:alnum:]]?[_]?)+ {
	yylval.str = strdup(yytext);
	return STR;
}

. { printf("Unexpected character %s@%d\n", yytext,yylineno); exit(0);}
%%

