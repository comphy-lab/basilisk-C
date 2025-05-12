/**
The problme is in [src/ast/Makefile](http://basilisk.fr/src/ast/Makefile#22)

*/
grammar: grammar.c
	cc -Wall grammar.c -o grammar
  
/**
My old `cc` binary fail to compile grammar.c. I think it should be
*/

grammar: grammar.c
	$(CC) -Wall grammar.c -o grammar
