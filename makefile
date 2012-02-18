zice: syntactic.tab.c lex.yy.c main.c solution.c csparse.c analysis.c sparse_analysis.c utility.c 
	gcc -g lex.yy.c syntactic.tab.c solution.c main.c csparse.c analysis.c sparse_analysis.c utility.c -lm -o zice

syntactic.tab.c: syntactic.y types.h
	bison -vt --defines=syntactic.tab.h syntactic.y

lex.yy.c: lexical.l syntactic.tab.h types.h
	flex -i lexical.l

clean:
	rm -f syntactic.tab.c lex.yy.c debug zice syntactic.output syntactic.tab.h

cloc: clean 
	cloc lexical.l syntactic.y solution.* types.h main.c analysis.* sparse_analysis.* utility.* 
