/**
# Lexer for Basilisk C

Closely based on the [C99 lexer](c.lex). */

%option bison-bridge noyywrap yylineno nounput

%e  1019
%p  2807
%n  371
%k  284
%a  1213
%o  1117

O   [0-7]
D   [0-9]
NZ  [1-9]
L   [a-zA-Z_]
A   [a-zA-Z_0-9]
H   [a-fA-F0-9]
HP  (0[xX])
E   ([Ee][+-]?{D}+)
P   ([Pp][+-]?{D}+)
FS  (f|F|l|L)
IS  (((u|U)(l|L|ll|LL)?)|((l|L|ll|LL)(u|U)?))
CP  (u|U|L)
SP  (u8|u|U|L)
ES  (\\(['"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
WS  [ \t\v\n\f]
STRING \"([^"\\\n]|{ES})*\"

%{
#include <stdio.h>
#include <assert.h>
#include "parser.h"
#include "basilisk.h"

#define YYSTYPE Ast *
#define YY_DECL int yylex (YYSTYPE * yylval_param, AstRoot * parse)
 
static void comment (void);
static void preproc (void);
static void bpreproc (void);
static void ompreproc (void);
static void file_line (AstRoot * parse, const char * text);
static int  check_type (AstRoot * parse);

static Ast * new_ast (AstRoot * parse,
		      int token, int line, char * start, char * end)
{
  Allocator * alloc = parse->alloc;
  Ast * n = allocate (alloc, sizeof(AstTerminal));
  memset (n, 0, sizeof(AstTerminal));
  n->sym = token_symbol (token);
  AstTerminal * t = ast_terminal (n);
  t->start = start;
  t->after = end;
  t->file = parse->file;
  while (start != end) {
    if (*start == '\n')
      line--;
    start++;
  }
  t->line = line;
  return n;
}

#define CAST()								\
  *yylval = new_ast (parse, yytext[0], yylineno, yytext, yytext);	\
  return yytext[0];
  
#define SAST(t) 							\
  *yylval = new_ast (parse, t, yylineno, yytext, yytext + strlen(yytext) - 1); \
  return t;
  
%}
	 
%%

"/*"                                    { comment(); }
"//".*                                  { /* consume //-comment */ }
^[ \t]*#[ \t]+[0-9]+[ \t]+{STRING}.*    { file_line (parse, yytext); }
^[ \t]*#	                        { preproc(); }
^[ \t]*@[ \t]*def[ \t].*                { bpreproc(); }
^[ \t]*@.*
^[ \t]*OMP[ \t]*\(	                { ompreproc(); }
	 
"auto"					{ SAST(AUTO); }
"break"					{ SAST(BREAK); }
"case"					{ SAST(CASE); }
"char"					{ SAST(CHAR); }
"const"					{ SAST(CONST); }
"continue"				{ SAST(CONTINUE); }
"default"				{ SAST(DEFAULT); }
"do"					{ SAST(DO); }
"double"				{ SAST(DOUBLE); }
"else"					{ SAST(ELSE); }
"enum"					{ SAST(ENUM); }
"extern"				{ SAST(EXTERN); }
"float"					{ SAST(FLOAT); }
"for"					{ SAST(FOR); }
"goto"					{ SAST(GOTO); }
"if"					{ SAST(IF); }
"inline"				{ SAST(INLINE); }
"int"					{ SAST(INT); }
"long"					{ SAST(LONG); }
"register"				{ SAST(REGISTER); }
"restrict"				{ SAST(RESTRICT); }
"return"				{ SAST(RETURN); }
"short"					{ SAST(SHORT); }
"signed"				{ SAST(SIGNED); }
"sizeof"				{ SAST(SIZEOF); }
"static"				{ SAST(STATIC); }
"struct"				{ SAST(STRUCT); }
"switch"				{ SAST(SWITCH); }
"typedef"				{ SAST(TYPEDEF); }
"union"					{ SAST(UNION); }
"unsigned"				{ SAST(UNSIGNED); }
"void"					{ SAST(VOID); }
"volatile"				{ SAST(VOLATILE); }
"while"					{ SAST(WHILE); }
"_Alignas"                              { SAST(ALIGNAS); }
"_Alignof"                              { SAST(ALIGNOF); }
"_Atomic"                               { SAST(ATOMIC); }
"_Bool"                                 { SAST(BOOL); }
"_Complex"                              { SAST(COMPLEX); }
"complex"                               { SAST(COMPLEX); }
"_Generic"                              { SAST(GENERIC); }
"_Imaginary"                            { SAST(IMAGINARY); }
"_Noreturn"                             { SAST(NORETURN); }
"_Static_assert"                        { SAST(STATIC_ASSERT); }
"_Thread_local"                         { SAST(THREAD_LOCAL); }
"__func__"                              { SAST(FUNC_NAME); }

                    /* Basilisk C tokens */

"new"{WS}+("vertex"{WS}+)?"scalar"      { SAST(NEW_FIELD); }
"new"{WS}+("face"{WS}+)?"vector"        { SAST(NEW_FIELD); }
"new"{WS}+("symmetric"{WS}+)?"tensor"   { SAST(NEW_FIELD); }
"vertex"{WS}+"scalar"                   { SAST(TYPEDEF_NAME); }
"face"{WS}+"vector"                     { SAST(TYPEDEF_NAME); }
"symmetric"{WS}+"tensor"                { SAST(TYPEDEF_NAME); }
"(const)"                               { SAST(MAYBECONST); }
"trace"			                { SAST(TRACE); }
"reduction"			        { SAST(REDUCTION); }

"foreach_blockf" |
"foreach_block" |
"foreach_child" |
"foreach_neighbor"                      { SAST(FOREACH_INNER); }

"foreach_dimension"			{ SAST(FOREACH_DIMENSION); }

"foreach" |
"foreach_"{L}{A}*                       { SAST(FOREACH); }

                    /* GCC extensions */
	   
"__attribute__"{WS}*\(                  { ompreproc(); }
	   
                    /* End of GCC extensions */

{L}{A}*					{ SAST(check_type (parse)); }

{HP}{H}+{IS}?				{ SAST(I_CONSTANT); }
{NZ}{D}*{IS}?				{ SAST(I_CONSTANT); }
"0"{O}*{IS}?				{ SAST(I_CONSTANT); }
{CP}?"'"([^'\\\n]|{ES})+"'"		{ SAST(I_CONSTANT); }

{D}+{E}{FS}?				{ SAST(F_CONSTANT); }
{D}*"."{D}+{E}?{FS}?			{ SAST(F_CONSTANT); }
{D}+"."{E}?{FS}?			{ SAST(F_CONSTANT); }
{HP}{H}+{P}{FS}?			{ SAST(F_CONSTANT); }
{HP}{H}*"."{H}+{P}{FS}?			{ SAST(F_CONSTANT); }
{HP}{H}+"."{P}{FS}?			{ SAST(F_CONSTANT); }

({SP}?\"([^"\\\n]|{ES})*\"{WS}*)+	{ SAST(STRING_LITERAL); }

"..."					{ SAST(ELLIPSIS); }
">>="					{ SAST(RIGHT_ASSIGN); }
"<<="					{ SAST(LEFT_ASSIGN); }
"+="					{ SAST(ADD_ASSIGN); }
"-="					{ SAST(SUB_ASSIGN); }
"*="					{ SAST(MUL_ASSIGN); }
"/="					{ SAST(DIV_ASSIGN); }
"%="					{ SAST(MOD_ASSIGN); }
"&="					{ SAST(AND_ASSIGN); }
"^="					{ SAST(XOR_ASSIGN); }
"|="					{ SAST(OR_ASSIGN); }
">>"					{ SAST(RIGHT_OP); }
"<<"					{ SAST(LEFT_OP); }
"++"					{ SAST(INC_OP); }
"--"					{ SAST(DEC_OP); }
"->"					{ SAST(PTR_OP); }
"&&"					{ SAST(AND_OP); }
"||"					{ SAST(OR_OP); }
"<="					{ SAST(LE_OP); }
">="					{ SAST(GE_OP); }
"=="					{ SAST(EQ_OP); }
"!="					{ SAST(NE_OP); }
";"					{ CAST(); }
("{"|"<%")				{ CAST(); }
("}"|"%>")				{ CAST(); }
","					{ CAST(); }
":"					{ CAST(); }
"="					{ CAST(); }
"("					{ CAST(); }
")"					{ CAST(); }
("["|"<:")				{ CAST(); }
("]"|":>")				{ CAST(); }
"."					{ CAST(); }
"&"					{ CAST(); }
"!"					{ CAST(); }
"~"					{ CAST(); }
"-"					{ CAST(); }
"+"					{ CAST(); }
"*"					{ CAST(); }
"/"					{ CAST(); }
"%"					{ CAST(); }
"<"					{ CAST(); }
">"					{ CAST(); }
"^"					{ CAST(); }
"|"					{ CAST(); }
"?"					{ CAST(); }

{WS}					{ /* whitespace separates tokens */ }
.					{ /* discard bad characters */ }
	 
%%

static void comment (void)
{
  int c;  
  while ((c = input()) != 0)
    if (c == '*') {
      while ((c = input()) == '*')
	;
      
      if (c == '/')
	return;
      
      if (c == 0)
	break;
    }
  //  yyerror ("unterminated comment");
}

static void preproc (void)
{
  int c, c1 = 0;
  while ((c = input()) != 0) {
    if (c == '\n' && c1 != '\\')      
      return;
    c1 = c;
  }
  //  yyerror ("unterminated preprocessor directive");
}

static void bpreproc (void)
{
  int c;
  while ((c = input()) != 0)
    if (c == '@')
      return;
  //  yyerror ("unterminated @def");
}

static void ompreproc (void)
{
  int c, scope = 1;
  while ((c = input()) != 0) {
    if (c == '(')
      scope++;
    else if (c == ')') {
      scope--;
      if (scope == 0)
	return;
    }
  }
  //  yyerror ("unterminated OMP");
}

static void file_line (AstRoot * parse, const char * text)
{
  char * s = strchr (text, '#') + 1;
  yylineno = atoi(s) - 1;
  s = strchr (s, '"') + 1;
  char * end = strchr (s, '"');
  parse->file = allocate (parse->alloc, end - s + 1);
  strncpy ((char *) parse->file, s, end - s);
  //  fprintf (stderr, "%s: \"%s\" %d\n", text, file, yylineno);
}

static int check_type (AstRoot * parse)
{
  if (parse->type_already_specified)
    return IDENTIFIER;
  
  Ast * declaration = ast_identifier_declaration (parse->stack, yytext);
  if (declaration) {
    if (ast_is_typedef (declaration))
      return TYPEDEF_NAME;
    return IDENTIFIER;
  }

  return IDENTIFIER;
}

void lexer_setup (char * buffer, size_t len)
{
  yylineno = 1;
  yy_scan_buffer (buffer, len);
}
