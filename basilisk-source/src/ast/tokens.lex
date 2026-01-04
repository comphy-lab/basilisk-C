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
#include "symbols.h"
	   
#define YYSTYPE Ast *
#define YY_DECL int yylex (YYSTYPE * yylval_param, AstRoot * parse)
 
static void comment (void);
static void preproc (void);
static void bpreproc (void);
static void ompreproc (void);
static void file_line (AstRoot * parse, const char * text);
static int  check_type (AstRoot * parse, bool call);

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
	 
                    /* Basilisk C tokens */

"new"{WS}+("vertex"{WS}+)?"scalar"      { SAST(NEW_FIELD); }
"new"{WS}+("face"{WS}+)?"vector"        { SAST(NEW_FIELD); }
"new"{WS}+("symmetric"{WS}+)?"tensor"   { SAST(NEW_FIELD); }
"vertex"{WS}+"scalar"                   { SAST(TYPEDEF_NAME); }
"face"{WS}+"vector"                     { SAST(TYPEDEF_NAME); }
"symmetric"{WS}+"tensor"                { SAST(TYPEDEF_NAME); }
"(const)"                               { SAST(MAYBECONST); }
"trace"			                { SAST(TRACE); }

                    /* GCC extensions */
	   
"__attribute__"{WS}*\(                  { ompreproc(); }
	   
                    /* End of GCC extensions */

{L}{A}*/{WS}*"("			{ SAST(check_type (parse, true)); }
{L}{A}*			                { SAST(check_type (parse, false)); }

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

"{...}"					{ SAST(ELLIPSIS_MACRO); }
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

static int check_type (AstRoot * parse, bool call)
{
  typedef struct {
    const char * name;
    int type;
  } Type;
  static Type types[] = {
    {"auto",	       AUTO},
    {"break",	       BREAK},
    {"case",	       CASE},
    {"char",	       CHAR},
    {"const",	       CONST},
    {"continue",       CONTINUE},
    {"default",	       DEFAULT},
    {"do",	       DO},
    {"double",	       DOUBLE},
    {"else",	       ELSE},
    {"enum",	       ENUM},
    {"extern",	       EXTERN},
    {"float",	       FLOAT},
    {"for",	       FOR},
    {"goto",	       GOTO},
    {"if",	       IF},
    {"inline",	       INLINE},
    {"int",	       INT},
    {"long",	       LONG},
    {"register",       REGISTER},
    {"restrict",       RESTRICT},
    {"return",	       RETURN},
    {"short",	       SHORT},
    {"signed",	       SIGNED},
    {"sizeof",	       SIZEOF},
    {"static",	       STATIC},
    {"struct",	       STRUCT},
    {"switch",	       SWITCH},
    {"typedef",	       TYPEDEF},
    {"union",	       UNION},
    {"unsigned",       UNSIGNED},
    {"void",	       VOID},
    {"volatile",       VOLATILE},
    {"while",          WHILE},
    {"_Alignas",       ALIGNAS},
    {"_Alignof",       ALIGNOF},
    {"_Atomic",        ATOMIC},
    {"_Bool",          BOOL},
    {"_Complex",       COMPLEX},
    {"complex",        COMPLEX},
    {"_Generic",       GENERIC},
    {"_Imaginary",     IMAGINARY},
    {"_Noreturn",      NORETURN},
    {"_Static_assert", STATIC_ASSERT},
    {"_Thread_local",  THREAD_LOCAL},
    {"__func__",       FUNC_NAME},

                        /* Basilisk C tokens */
    
    {"reduction", REDUCTION},
    {"foreach_dimension", FOREACH_DIMENSION},
    
    {NULL}};
  for (Type * t = types; t->name; t++)
    if (!strcmp (yytext, t->name))
      return t->type;

  if (!strncmp (yytext, "macro", 5)) {
    char * s;
    for (s = yytext + 5; *s != '\0' && *s >= '0' && *s <= '9'; s++);
    if (*s == '\0')
      return MACRODEF;
  }
  
  if (parse->type_already_specified)
    return IDENTIFIER;
  
  return ast_identifier_parse_type (parse->stack, yytext, call, parse->file, yylineno);
}

void lexer_setup (char * buffer, size_t len)
{
  yylineno = 1;
  yy_scan_buffer (buffer, len);
}
