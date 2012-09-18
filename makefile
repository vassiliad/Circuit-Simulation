zice: parser.o lex.o main.o hash_table.o components.o mna.o utility.o plot.o algebra.o transient.o csparse.o
	gcc -Wall -g *.o -lm  -o zice

csparse.o: csparse.h csparse.c
	gcc -Wall -g -c csparse.c -o csparse.o

parser.o: parser.tab.c
	gcc -Wall -g -c parser.tab.c -o parser.o

lex.o: lex.yy.c
	gcc -Wall -g -c lex.yy.c -o lex.o

main.o: main.c
	gcc -Wall -g -c main.c -o main.o

hash_table.o: hash_table.c
	gcc -Wall -g -c hash_table.c -o hash_table.o

components.o: components.c
	gcc -Wall -g -c components.c -o components.o

mna.o: mna.c mna.h
	gcc -Wall -g -c mna.c -o mna.o

utility.o: utility.c utility.h
	gcc -Wall -g -c utility.c -o utility.o

algebra.o: algebra.c algebra.h
	gcc -Wall -g -c algebra.c -o algebra.o

plot.o: plot.c plot.h
	gcc -Wall -g -c plot.c -o plot.o

transient.o: transient.c transient.h
	gcc -Wall -g -c transient.c -o transient.o

parser.tab.c: parser.y
	bison -vt --defines=parser.h parser.y

lex.yy.c: lexical.l parser.h
	flex -i lexical.l

clean:
	rm -f parser.tab.c lex.yy.c debug zice parser.output parser.h *.o

cloc: clean 
	cloc lexical.l parser.y  main.c components.* hash_table.* options.h utility.* mna.* solution.* transient.* algebra.*
