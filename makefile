zice: syntactic.tab.c lex.yy.c main.c solution.c
	gcc -g lex.yy.c syntactic.tab.c solution.c main.c -lm -o zice

debug: syntactic.tab.c lex.yy.c main.c
	gcc -Wall -lm -g -o debug lex.yy.c syntactic.tab.c main.c solution.c

syntactic.tab.c: syntactic.y types.h
	bison -vt --defines=syntactic.tab.h syntactic.y

lex.yy.c: lexical.l syntactic.tab.h types.h
	flex -i lexical.l

clean:
	rm -f syntactic.tab.c lex.yy.c debug zice syntactic.output syntactic.tab.h

cloc: clean
	cloc --exclude-lang=make,"Bourne Shell" ./