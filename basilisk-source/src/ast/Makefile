CFLAGS += -I.. -DBASILISK=\"$(subst /ast,,$(CURDIR))\"

OBJECTS = ast.o tokens.o basilisk.o translate.o allocator.o \
	faststack.o stencil.o types.o references.o kernels.o check.o interpreter/interpreter.o

TOPTARGETS = all clean check

SUBDIRS = interpreter

.PHONY: subdirs $(SUBDIRS) $(TOPTARGETS)

all: subdirs libast.a

subdirs: $(SUBDIRS)

$(TOPTARGETS): $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@ $(filter-out $(SUBDIRS),$(MAKECMDGOALS))

libast.a: $(OBJECTS)
	ar cr $@ $(OBJECTS)

# Uncomment the recipe below if you need to re-generate basilisk.c

# basilisk.c: basilisk.yacc
# 	bison -Dparse.trace -d -o basilisk.c basilisk.yacc
# 	sed -i 's/^ *YY_REDUCE_PRINT *(yyn);/ YY_REDUCE_PRINT (yyn); DEFAULT_ACTION (yyn);/' basilisk.c

symbols.h: basilisk.c Makefile
	echo "// Automatically generated from basilisk.yacc" > symbols.h
	@echo 'const char * symbol_name (int sym);' >> symbols.h
	@echo 'int token_symbol (int token);' >> symbols.h
	@awk '/^enum yysymbol_kind_t/{a=1}{if(a)print $$0;}/};/{a=0}' \
		basilisk.c | sed 's/YYSYMBOL_/sym_/g' >> symbols.h

grammar.h: basilisk.yacc grammar
	./grammar < basilisk.yacc  > grammar.h

grammar: grammar.c
	$(CC) $(CFLAGS) grammar.c -o grammar

# Uncomment the recipe below if you need to re-generate tokens.c

# tokens.c: tokens.lex
#	flex -o tokens.c tokens.lex

basilisk.o: ast.h parser.h

ast.o: ast.h parser.h symbols.h

stack.o: stack.h

faststack.o: stack.h symbols.h ast.h

tokens.o: ast.h basilisk.c

translate.o: ast.h symbols.h einstein_sum.h macro.h optional.h

stencil.o: ast.h symbols.h

allocator.o: allocator.h

check.o: ast.h symbols.h grammar.h

expr.o: ast.h symbols.h

types.o: ast.h symbols.h

references.o: ast.h symbols.h

kernels.o: ast.h symbols.h

clean:
	rm -f *.o *.a

qcc: qcc.c libast.a
	$(CC) $(CFLAGS) -o $@ $^

expr: expr.c ast.h libast.a
	$(CC) $(CFLAGS) -o $@ $< -L. -last -lm

