%option noyywrap
%option yylineno
%{
#include <stdio.h>
#include <ctype.h>
#include "components.h"
#include "parser.h"

#define zprint(...) 
#define mprint printf

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

".OPTIONS" return OPTIONS;
".OPTION"  return OPTIONS;

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
	yylval.integer = strtoul(yytext, 0, 10);
  zprint("Number_i: %d\n", yylval.integer);
	return INTEGER;
}

[-+]?[[:digit:]]+(\.[[:digit:]]+(e[-+]?[[:digit:]]+)?|(\.[[:digit:]]+)?(e[-+]?[[:digit:]]+)) {
	yylval.dbl = strtod(yytext,NULL);
  zprint("Number_d: %g\n", yylval.dbl);
	return DOUBLE;
}



V[ \t]*\([ \t]*[[:digit:]]+[ \t]*\) {
  char * str;
  int i,j;

  j = 0;
  str = strchr(yytext, '(')+1;
  char *s = strdup(str);

  j = 0;
  for (i=0; str[i] != ')'; i++ ) {
    if ( str[i] == ' ' || str[i] == '\t' )
      continue;
    s[j++] = str[i];
  }

  s[j] = 0;
  yylval.string = s;
	return PLOT_V;
}

\n+ {
	return NEW_LINE;
}

[\ \r\t]+ { }

[_]?[[:alpha:]]+([[:alnum:]]?[_]?)+ {
#if 0
  int i = 0 ;
  int j = 0;

  char *s = strdup(yytext);

  for (i=0; i<strlen(yytext); i++ ) {
    if ( yytext[i] == ' ' || yytext[i] == '\t' )
      continue;
    s[j++] = yytext[i];

  }
  s[j] = 0;
  yylval.string = s;
  zprint("str: %s\n", s);
  mprint("%p %p\n", yytext, yylval.string);
#endif
  yylval.string = strdup(yytext);
	return STRING;
}

. { printf("Unexpected character %s@%d\n", yytext,yylineno); exit(0);}
%%

