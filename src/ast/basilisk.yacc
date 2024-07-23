/**
# Yacc/Bison Grammar for Basilisk C

The grammar is closely based on the [C99 grammar](c.yacc).

The main difference with the C99 grammar is the addition of the
Basilisk C syntax and the further subdivision of some the C99 entries,
to simplify parsing.

A classical difficulty with the C grammar is that it is "almost
context-free" but not quite. The problem comes from the different uses
in the grammar of identifiers and typedefs. This problem has been known
for a long time and complicates the design of a Yacc grammar for C
(see [Jourdan et al., 2017](#jourdan2017) for a comprehensive review
and solutions).

In the following we use a ["lexical
tie-in"](https://www.gnu.org/software/bison/manual/html_node/Lexical-Tie_002dins.html)
as described in the [bison
manual](https://www.gnu.org/software/bison/manual/html_node/index.html)
to "solve" the issue. This corresponds with the `type_not_specified`
entry in the grammar. I put quotes around "solve" because solving this
context issue truly consistently is very difficult, as described in
[Jourdan et al., 2017](#jourdan2017). The level of robustness of this
approach should be appropriate for the vast majority of use cases of
Basilisk, however.

The C and Basilisk C tokens are in [tokens.lex]().

## References

~~~bib
@hal{jourdan2017, hal-01633123}

@article{scarpazza2007,
author = {Scarpazza, D.P.},
year = {2007},
pages = {48-55},
title = {Practical parsing for ANSI C},
volume = {32},
journal = {Dr. Dobb's Journal}
}

Jim Roskind's C grammar
https://blog.robertelder.org/jim-roskind-grammar/c%2B%2Bgrammar2.0.tar.Z
~~~
*/

%param { AstRoot * parse }
%parse-param { Ast * root }
%define api.pure full
%define api.value.type { Ast * }

%{
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include "parser.h"

static Ast * ast_reduce (Allocator * alloc, int sym, Ast ** children, int n);
#define DEFAULT_ACTION(yyn)					\
  yyval = ast_reduce ((Allocator *)parse->alloc, yyr1[yyn], yyvsp, yyr2[yyn])
static int yyparse (AstRoot * parse, Ast * root);
%}

/**
## C99 tokens */

%token	IDENTIFIER I_CONSTANT F_CONSTANT STRING_LITERAL FUNC_NAME SIZEOF
%token	PTR_OP INC_OP DEC_OP LEFT_OP RIGHT_OP LE_OP GE_OP EQ_OP NE_OP
%token	AND_OP OR_OP MUL_ASSIGN DIV_ASSIGN MOD_ASSIGN ADD_ASSIGN
%token	SUB_ASSIGN LEFT_ASSIGN RIGHT_ASSIGN AND_ASSIGN
%token	XOR_ASSIGN OR_ASSIGN
%token	TYPEDEF_NAME ENUMERATION_CONSTANT

%token	TYPEDEF EXTERN STATIC AUTO REGISTER INLINE
%token	CONST RESTRICT VOLATILE
%token	BOOL CHAR SHORT INT LONG SIGNED UNSIGNED FLOAT DOUBLE VOID
%token	COMPLEX IMAGINARY 
%token	STRUCT UNION ENUM ELLIPSIS

%token	CASE DEFAULT IF ELSE SWITCH WHILE DO FOR GOTO CONTINUE BREAK RETURN

%token	ALIGNAS ALIGNOF ATOMIC GENERIC NORETURN STATIC_ASSERT THREAD_LOCAL

/**
## Basilisk C tokens */

%token  MAYBECONST NEW_FIELD TRACE
%token  FOREACH FOREACH_INNER FOREACH_DIMENSION
%token  REDUCTION

/**
## Grammar

Note that 9 shift/reduce conflicts and 5 reduce/reduce conflicts are
expected. */

%start root

%%

translation_unit
        : external_declaration
	| translation_unit external_declaration	
        | translation_unit error ';'  { $2->sym = YYSYMBOL_YYerror; }
        | translation_unit error '}'  { $2->sym = YYSYMBOL_YYerror; }
        | translation_unit error ')'  { $2->sym = YYSYMBOL_YYerror; }
        ;

primary_expression
        : IDENTIFIER
	| constant
	| string
	| '(' expression_error ')'
	| generic_selection
	;

expression_error
        : expression
	| error         { $1->sym = YYSYMBOL_YYerror; }
	;

constant
	: I_CONSTANT		/* includes character_constant */
	| F_CONSTANT
	| ENUMERATION_CONSTANT	/* after it has been defined as such */
	;

enumeration_constant		/* before it has been defined as such */
	: IDENTIFIER
	;

string
	: STRING_LITERAL
	| FUNC_NAME
	;

generic_selection
        : GENERIC '(' assignment_expression ',' generic_assoc_list ')'
	;

generic_assoc_list
	: generic_association
	| generic_assoc_list ',' generic_association
	;

generic_association
	: type_name ':' assignment_expression
	| DEFAULT ':' assignment_expression
	;

postfix_expression
	: primary_expression
	| function_call
	| array_access
        | postfix_expression '.' member_identifier
        | postfix_expression PTR_OP member_identifier
	| postfix_expression INC_OP
	| postfix_expression DEC_OP
	| '(' type_name ')' postfix_initializer
	;

postfix_initializer
        : '{' initializer_list '}'
	| '{' initializer_list ',' '}'
	;
	
array_access
        : postfix_expression '[' ']' /* Basilisk C extension */
	| postfix_expression '[' '*' ']' /* Basilisk C extension */
        | postfix_expression '[' expression ']'
	;
	
function_call
        : postfix_expression '(' ')'
	| postfix_expression '(' argument_expression_list ')'
        ;

member_identifier
        : generic_identifier
	;

argument_expression_list
	: argument_expression_list_item
	| argument_expression_list ',' argument_expression_list_item
	| argument_expression_list ',' /* Basilisk C extension */
	;

argument_expression_list_item
        : assignment_expression
	| postfix_initializer /* Basilisk C extension */
	;

unary_expression
	: postfix_expression
	| INC_OP unary_expression
	| DEC_OP unary_expression
	| unary_operator cast_expression
	| SIZEOF unary_expression
	| SIZEOF '(' type_name ')'
	| ALIGNOF '(' type_name ')'
	| new_field /* Basilisk C extension */
	;

unary_operator
	: '&'
	| '*'
	| '+'
	| '-'
	| '~'
	| '!'
	;

cast_expression
	: unary_expression
	| '(' type_name ')' cast_expression
	;

/* fixme: simplify using precedence rules */

multiplicative_expression
	: cast_expression
	| multiplicative_expression '*' cast_expression
	| multiplicative_expression '/' cast_expression
	| multiplicative_expression '%' cast_expression
	;

additive_expression
	: multiplicative_expression
	| additive_expression '+' multiplicative_expression
	| additive_expression '-' multiplicative_expression
	;

shift_expression
	: additive_expression
	| shift_expression LEFT_OP additive_expression
	| shift_expression RIGHT_OP additive_expression
	;

relational_expression
	: shift_expression
	| relational_expression '<' shift_expression
	| relational_expression '>' shift_expression
	| relational_expression LE_OP shift_expression
	| relational_expression GE_OP shift_expression
	;

equality_expression
	: relational_expression
	| equality_expression EQ_OP relational_expression
	| equality_expression NE_OP relational_expression
	;

and_expression
	: equality_expression
	| and_expression '&' equality_expression
	;

exclusive_or_expression
	: and_expression
	| exclusive_or_expression '^' and_expression
	;

inclusive_or_expression
	: exclusive_or_expression
	| inclusive_or_expression '|' exclusive_or_expression
	;

logical_and_expression
	: inclusive_or_expression
	| logical_and_expression AND_OP inclusive_or_expression
	;

logical_or_expression
	: logical_and_expression
	| logical_or_expression OR_OP logical_and_expression
	;

conditional_expression
	: logical_or_expression
	| logical_or_expression '?' expression ':' conditional_expression
	;

assignment_expression
	: conditional_expression
	| unary_expression assignment_operator assignment_expression
	| unary_expression assignment_operator postfix_initializer /* Basilisk C extension */
	| TYPEDEF_NAME assignment_operator assignment_expression /* Basilisk C extension */
	| TYPEDEF_NAME assignment_operator postfix_initializer /* Basilisk C extension */
	;

assignment_operator
	: '='
	| MUL_ASSIGN
	| DIV_ASSIGN
	| MOD_ASSIGN
	| ADD_ASSIGN
	| SUB_ASSIGN
	| LEFT_ASSIGN
	| RIGHT_ASSIGN
	| AND_ASSIGN
	| XOR_ASSIGN
	| OR_ASSIGN
	;

expression
	: assignment_expression
	| expression ',' assignment_expression
	;

constant_expression
	: conditional_expression	/* with constraints */
	;

declaration
        : declaration_specifiers ';' type_not_specified                        { ast_push_declaration (parse->stack, $$); }
	| declaration_specifiers init_declarator_list ';' type_not_specified   { ast_push_declaration (parse->stack, $$); }
	| static_assert_declaration
	;

declaration_specifiers
        : storage_class_specifier declaration_specifiers
	| storage_class_specifier
	| type_specifier declaration_specifiers
	| type_specifier
	| type_qualifier declaration_specifiers
	| type_qualifier
	| function_specifier declaration_specifiers
	| function_specifier
	| alignment_specifier declaration_specifiers
	| alignment_specifier
	;

init_declarator_list
	: init_declarator
	| init_declarator_list ',' init_declarator
	;

init_declarator
	: declarator '=' initializer
	| declarator
	;

storage_class_specifier
	: TYPEDEF	/* identifiers must be flagged as TYPEDEF_NAME */
	| EXTERN
	| STATIC
	| THREAD_LOCAL
	| AUTO
	| REGISTER
	| TRACE /* Basilisk C extension */
	;

type_specifier
        : types    { parse->type_already_specified = true; }
        ;

types
        : VOID
	| CHAR
	| SHORT
	| INT
	| LONG
	| FLOAT
	| DOUBLE
	| SIGNED
	| UNSIGNED
	| BOOL
	| COMPLEX
	| IMAGINARY /* non-mandated extension */
	| atomic_type_specifier
	| struct_or_union_specifier
	| enum_specifier
	| TYPEDEF_NAME		/* after it has been defined as such */
	;

struct_or_union_specifier
	: struct_or_union '{' struct_declaration_list '}'
	| struct_or_union generic_identifier '{' struct_declaration_list '}'
	| struct_or_union generic_identifier
	| struct_or_union '{' error '}'                     { $3->sym = YYSYMBOL_YYerror; }
	| struct_or_union generic_identifier '{' error '}'  { $4->sym = YYSYMBOL_YYerror; }
	;

struct_or_union
	: STRUCT
	| UNION
	;

struct_declaration_list
	: struct_declaration
	| struct_declaration_list struct_declaration
	;

struct_declaration
        : specifier_qualifier_list ';' type_not_specified      /* for anonymous struct/union */
	| specifier_qualifier_list struct_declarator_list ';' type_not_specified
	| static_assert_declaration
	;

specifier_qualifier_list
	: type_specifier specifier_qualifier_list
	| type_specifier
	| type_qualifier specifier_qualifier_list
	| type_qualifier
	;

struct_declarator_list
	: struct_declarator
	| struct_declarator_list ',' struct_declarator
	;

struct_declarator
	: ':' constant_expression
	| declarator ':' constant_expression
	| declarator
	;

enum_specifier
	: ENUM '{' enumerator_list '}'
	| ENUM '{' enumerator_list ',' '}'
	| ENUM generic_identifier '{' enumerator_list '}'
	| ENUM generic_identifier '{' enumerator_list ',' '}'
	| ENUM generic_identifier
	;

enumerator_list
	: enumerator
	| enumerator_list ',' enumerator
	;

enumerator	/* identifiers must be flagged as ENUMERATION_CONSTANT */
	: enumeration_constant '=' constant_expression
	| enumeration_constant
	;

atomic_type_specifier
	: ATOMIC '(' type_name ')'
	;

type_qualifier
	: CONST
	| RESTRICT
	| VOLATILE
	| ATOMIC
	| MAYBECONST
	;

function_specifier
	: INLINE
	| NORETURN
	;

alignment_specifier
	: ALIGNAS '(' type_name ')'
	| ALIGNAS '(' constant_expression ')'
	;

declarator
        : pointer direct_declarator type_not_specified
	| direct_declarator type_not_specified
	;

direct_declarator
        : generic_identifier
	| '(' declarator ')'
	| direct_declarator '[' ']'
	| direct_declarator '[' '*' ']'
	| direct_declarator '[' STATIC type_qualifier_list assignment_expression ']'
	| direct_declarator '[' STATIC assignment_expression ']'
	| direct_declarator '[' type_qualifier_list '*' ']'
	| direct_declarator '[' type_qualifier_list STATIC assignment_expression ']'
	| direct_declarator '[' type_qualifier_list assignment_expression ']'
	| direct_declarator '[' type_qualifier_list ']'
	| direct_declarator '[' assignment_expression ']'
	| direct_declarator '(' type_not_specified parameter_type_list ')'
	| direct_declarator '(' type_not_specified error ')'                               { $4->sym = YYSYMBOL_YYerror; }
	| direct_declarator '(' type_not_specified ')'
	| direct_declarator '(' type_not_specified identifier_list ')'
	;

generic_identifier
        : IDENTIFIER
	| TYPEDEF_NAME
	;

pointer
	: '*' type_qualifier_list pointer
	| '*' type_qualifier_list
	| '*' pointer
	| '*'
	;

type_qualifier_list
	: type_qualifier
	| type_qualifier_list type_qualifier
	;

parameter_type_list
	: parameter_list ',' ELLIPSIS
	| parameter_list
	;

parameter_list
        : parameter_declaration type_not_specified
	| parameter_list ',' parameter_declaration type_not_specified
	;

parameter_declaration
	: declaration_specifiers declarator
	| declaration_specifiers declarator '=' initializer /* Basilisk C extension */
	| declaration_specifiers abstract_declarator
	| declaration_specifiers
	;

identifier_list
	: IDENTIFIER
	| identifier_list ',' IDENTIFIER
	;

type_name
        : specifier_qualifier_list abstract_declarator type_not_specified
	| specifier_qualifier_list type_not_specified
	;

abstract_declarator
	: pointer direct_abstract_declarator
	| pointer
	| direct_abstract_declarator
	;

direct_abstract_declarator
	: '(' abstract_declarator ')'
	| '[' ']'
	| '[' '*' ']'
	| '[' STATIC type_qualifier_list assignment_expression ']'
	| '[' STATIC assignment_expression ']'
	| '[' type_qualifier_list STATIC assignment_expression ']'
	| '[' type_qualifier_list assignment_expression ']'
	| '[' type_qualifier_list ']'
	| '[' assignment_expression ']'
	| direct_abstract_declarator '[' ']'
	| direct_abstract_declarator '[' '*' ']'
	| direct_abstract_declarator '[' STATIC type_qualifier_list assignment_expression ']'
	| direct_abstract_declarator '[' STATIC assignment_expression ']'
	| direct_abstract_declarator '[' type_qualifier_list assignment_expression ']'
	| direct_abstract_declarator '[' type_qualifier_list STATIC assignment_expression ']'
	| direct_abstract_declarator '[' type_qualifier_list ']'
	| direct_abstract_declarator '[' assignment_expression ']'
	| '(' ')'
	| '(' parameter_type_list ')'
	| direct_abstract_declarator '(' type_not_specified ')'
	| direct_abstract_declarator '(' type_not_specified parameter_type_list ')'
	;

type_not_specified
        :  { parse->type_already_specified = false; }
        ;

initializer
	: '{' initializer_list '}'
	| '{' initializer_list ',' '}'
	| assignment_expression
	;

initializer_list
	: designation initializer
	| initializer
	| initializer_list ',' designation initializer
	| initializer_list ',' initializer
	;

designation
	: designator_list '='
	;

designator_list
	: designator
	| designator_list designator
	;

designator
	: '[' constant_expression ']'
	| '.' generic_identifier
	;

static_assert_declaration
	: STATIC_ASSERT '(' constant_expression ',' STRING_LITERAL ')' ';'
	;

statement
	: labeled_statement
	| compound_statement
	| expression_statement
	| selection_statement
	| iteration_statement
	| jump_statement
        | basilisk_statements /* Basilisk C extension */
	| error ';'  { $1->sym = YYSYMBOL_YYerror; }
	;

labeled_statement
	: generic_identifier ':' statement
	| CASE constant_expression ':' statement
	| DEFAULT ':' statement
	;

compound_statement
	: '{' '}'
	|
	'{'                { stack_push (parse->stack, &($1)); $$->sym = YYSYMBOL_YYUNDEF; }
	block_item_list
	'}'	           { ast_pop_scope (parse->stack, $1); }
	;

block_item_list
	: block_item
	| block_item_list block_item
	;

block_item
	: declaration
	| statement
	;

expression_statement
	: ';'
	| expression ';'
	;

selection_statement
        : IF '(' expression_error ')' statement ELSE statement
        | IF '(' expression_error ')' statement
	| SWITCH '(' expression_error ')' statement
	;

for_scope
        : FOR { stack_push (parse->stack, &($$)); }
        ;

iteration_statement
        : WHILE '(' expression ')' statement                                            
	| DO statement WHILE '(' expression ')' ';'
	| for_scope '(' expression_statement expression_statement ')' statement
	            { ast_pop_scope (parse->stack, $1); }
	| for_scope '(' expression_statement expression_statement expression ')' statement
		    { ast_pop_scope (parse->stack, $1); }
	| for_declaration_statement
	;

for_declaration_statement
        : for_scope '(' declaration expression_statement ')' statement
	            { ast_pop_scope (parse->stack, $1); }	
	| for_scope '(' declaration expression_statement expression ')' statement
	            { ast_pop_scope (parse->stack, $1); }	
	;

jump_statement
        : GOTO generic_identifier ';'
	| CONTINUE ';'
	| BREAK ';'
	| RETURN ';'
	| RETURN expression ';'
	;

external_declaration
	: function_definition
	| declaration
	| macro_statement /* Basilisk C extension */
	| event_definition /* Basilisk C extension */
	| boundary_definition /* Basilisk C extension */
	| external_foreach_dimension /* Basilisk C extension */
	| attribute /* Basilisk C extension */
	| error compound_statement              { $1->sym = YYSYMBOL_YYerror; }
	;

function_declaration
        : declaration_specifiers declarator { ast_push_function_definition (parse->stack, $2);  }
	;
	
function_definition
        : function_declaration declaration_list compound_statement
                              	                { ast_pop_scope (parse->stack, $1->child[1]); }
	| function_declaration compound_statement
	                                        { ast_pop_scope (parse->stack, $1->child[1]); }
	;

declaration_list
	: declaration
	| declaration_list declaration
	;

/**
## Basilisk C grammar extensions */

basilisk_statements
        : macro_statement
        | foreach_statement
	| foreach_inner_statement
	| foreach_dimension_statement
	| forin_declaration_statement
	| forin_statement
	;

macro_statement
        : function_call compound_statement
        ;

foreach_statement
        : FOREACH '(' ')' statement
	| FOREACH '(' foreach_parameters ')' statement
	;

foreach_parameters
        : foreach_parameter
	| foreach_parameters ',' foreach_parameter
	;

foreach_parameter
        : assignment_expression
        | reduction_list
	;

reduction_list
        : reduction
	| reduction_list reduction
	;

reduction
        : REDUCTION '(' reduction_operator  ':' reduction_array ')'
	;

reduction_operator
        : generic_identifier
	| '+'
	| OR_OP
	;

reduction_array
        : generic_identifier
	| generic_identifier '[' ':' expression ']'
	;

foreach_inner_statement
        : FOREACH_INNER '(' ')' statement
	| FOREACH_INNER '(' expression ')' statement
	;

foreach_dimension_statement
        : FOREACH_DIMENSION '(' ')' statement
        | FOREACH_DIMENSION '(' I_CONSTANT ')' statement
	;

forin_declaration_statement
        : for_scope '(' declaration_specifiers declarator IDENTIFIER forin_arguments ')' statement
	            { ast_pop_scope (parse->stack, $1); }
        ;

forin_statement
        : for_scope '(' expression IDENTIFIER forin_arguments ')' statement
	            { ast_pop_scope (parse->stack, $1); }
	;

forin_arguments
        : expression
	| postfix_initializer
	;

event_definition
        : generic_identifier generic_identifier '(' event_parameters ')' statement
	;

event_parameters
        : event_parameter
	| event_parameters ',' event_parameter
	| event_parameters ';' event_parameter
        ;

event_parameter
        : conditional_expression
        | unary_expression assignment_operator conditional_expression
        | unary_expression assignment_operator postfix_initializer
        ;

boundary_definition
        : assignment_expression ';'
	;

external_foreach_dimension
        : FOREACH_DIMENSION '(' ')' function_definition
        | FOREACH_DIMENSION '(' I_CONSTANT ')' function_definition
	;

attribute
        : generic_identifier '{' struct_declaration_list '}'
	;

new_field
        : NEW_FIELD
	| NEW_FIELD '[' postfix_expression ']'
	;

root
        : translation_unit {
	  $$ = root;
	  $$->sym = yyr1[yyn];
	  $$->child = allocate ((Allocator *)parse->alloc, 2*sizeof(Ast *));
	  $$->child[0] = $1;
	  $1->parent = $$;
	  $$->child[1] = NULL;
        }
        ;

%%

/**
# Parsing functions */

/* Called by yyparse on error.  */
void
yyerror (AstRoot * parse, Ast * root, char const *s)
{
#if 0
  fprintf (stderr, "%d: %s near '", *line, s);
  char * s1 = *input - 1;
  while (!strchr("}{;\n", *s1)) s1--;
  s1++;
  while (strchr(" \t", *s1)) s1++;
  for (; s1 < *input; s1++)
    fputc (*s1, stderr);
  fputs ("'\n", stderr);
#endif
}

static char * copy_range (const char * start, const char * end, long offset)
{
  char * c = NULL;
  int len = end - start;
  if (len > 0) {
    char * s = c = malloc (len + 1);
    for (const char * i = start; i < end; i++, s++)
      *s = *(i + offset);
    *s = '\0';
  }
  return c;
}

static const char * copy_strings (const char * i, Ast * n, long offset)
{
  AstTerminal * t = ast_terminal (n);
  if (t) {
    t->before = copy_range (i, t->start, offset);
    if (t->start > i)
      i = t->start;

    t->start = copy_range (i, t->after + 1, offset);
    if (t->after + 1 > i)
      i = t->after + 1;
    t->after = NULL;
  }
  else
    for (Ast ** c = n->child; *c; c++)
      i = copy_strings (i, *c, offset);
  return i;
}

static void remove_child (Ast * c)
{
  if (!c->parent)
    return;
  Ast ** i = c->parent->child;
  for (; *i && *i != c; i++);
  assert (*i == c);
  for (; *i; i++)
    *i = *(i + 1);
  c->parent = NULL;
}

static Ast * ast_reduce (Allocator * alloc, int sym, Ast ** children, int n)
{
  Ast * ast = allocate (alloc, sizeof(Ast));
  memset (ast, 0, sizeof(Ast));
  ast->sym = sym;
  int ndef = 0;
  for (int i = 0; i < n; i++) {
    Ast * c = children[i + 1 - n];
    if (c->sym != YYSYMBOL_YYUNDEF && c->sym != YYSYMBOL_type_not_specified)
      ndef++;
  }
  ast->child = allocate (alloc, (ndef + 1)*sizeof(Ast *));
  ndef = 0;
  for (int i = 0; i < n; i++) {
    Ast * c = children[i + 1 - n];
    if (c->sym != YYSYMBOL_YYUNDEF && c->sym != YYSYMBOL_type_not_specified) {
      if (c->parent)
	remove_child (c);
      c->parent = ast;
      ast->child[ndef++] = c;
    }
    else
      assert (!c->parent);
  }
  ast->child[ndef] = NULL;
  return ast;
}

static Stack * stack_internalize (Stack * stack)
{
  Ast ** n;
  for (int i = 0; (n = stack_index (stack, i)); i++)
    if ((*n)->sym == YYSYMBOL_IDENTIFIER) {
      AstTerminal * t = ast_terminal (*n);
      char * after = t->start + strlen (t->start) - 1;
      if (t->after != NULL && t->after != after) {
	fprintf (stderr, "%s:%d: %s after: %s\n",
		 t->file, t->line, t->start, t->after);
	abort();
      }
      t->after = after;
    }
  return stack;
}

static void stack_externalize (Stack * stack)
{
  Ast ** n;
  for (int i = 0; (n = stack_index (stack, i)); i++)
    if ((*n)->sym == YYSYMBOL_IDENTIFIER) {
      AstTerminal * t = ast_terminal(*n);
      if (t->after != NULL) {
	if (t->after[1] != '\0') {
	  
	  /**
	  This is a declaration which has not been through
	  copy_strings() i.e. which is not connected to the root, due
	  to a syntax error which lead to the corresponding branch being
	  discarded. We set the symbol to UNDEF. */

	  t->start = t->before = NULL;
	  (*n)->sym = YYSYMBOL_YYUNDEF;
	}
	t->after = NULL;
      }
    }
}

AstRoot * ast_parse (const char * code, AstRoot * parent)
{
  AstRoot parse;
  parse.alloc = parent ? parent->alloc : new_allocator();
  parse.stack = parent ? parent->stack : stack_new (sizeof (Ast *));
  parse.file = allocate (parse.alloc, strlen ("<basilisk>") + 1);
  strcpy ((char *) parse.file, "<basilisk>");
  stack_internalize (parse.stack);
  parse.type_already_specified = false;
  extern void lexer_setup (char * buffer, size_t len);
  size_t len = strlen (code) + 1;
  char * buffer = malloc (len + 1);
  memcpy (buffer, code, len);
  buffer[len] = '\0';
  lexer_setup (buffer, len + 1);
  //  yydebug = 1;
  AstRoot * root = allocate ((Allocator *)parse.alloc, sizeof(AstRoot));
  memset (root, 0, sizeof(AstRoot));
  stack_push (parse.stack, &root);
  yyparse (&parse, (Ast *) root);
  if (((Ast *)root)->child) {
    const char * i = copy_strings (buffer, (Ast *) root, code - buffer);
    const char * end = i; while (*end != '\0') end++;
    root->after = copy_range (i, end, code - buffer);
    root->alloc = parent ? NULL : parse.alloc;
    root->stack = parent ? NULL : parse.stack;
    stack_externalize (parse.stack);
  }
  else {
    root = NULL;
    if (parent)
      stack_externalize (parse.stack);
    else {
      free_allocator (parse.alloc);
      stack_destroy (parse.stack);
    }
  }
  free (buffer);
  yylex_destroy();
  return root;
}

int token_symbol (int token)
{
  return YYTRANSLATE (token);
}

const char * symbol_name (int sym)
{
  return yytname[sym];
}
