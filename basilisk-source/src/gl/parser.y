%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
  
typedef struct _Allocator Allocator;  
  
struct _Allocator {
  void * m;
  long len, maxlen;
  Allocator * next;
};

static Allocator * new_allocator()
{
  Allocator * a = calloc (1, sizeof (Allocator));
  a->maxlen = 1 << 16;
  a->m = calloc (a->maxlen, 1);
  return a;
}

static void * allocate (Allocator * a, long size)
{
  Allocator * last = a;
  while (a && a->len + size >= a->maxlen)
    last = a, a = a->next;
  if (a == NULL)
    a = last->next = new_allocator();
  assert (a->len + size < a->maxlen);
  void * p = (void *)(((char *)a->m) + a->len);
  a->len += size;
  return p;
}

static void free_allocator (Allocator * a)
{
  free (a->m);
  if (a->next) free_allocator (a->next);
  free (a);
}
  
#include "parser.h"

static Node * new_value (Allocator * alloc, double val)
{
  Node * n = allocate (alloc, sizeof (Node));
  n->type = '1';
  n->d.value = val;
  return n;
}

static Node * new_expr (Allocator * alloc,
			char type, Node * e1, Node * e2, Node * e3)
{
  Node * n = allocate (alloc, sizeof (Node));
  n->type = type;
  n->e[0] = e1;
  n->e[1] = e2;
  n->e[2] = e3;
  return n;
}

static int yyparse (char ** input, Allocator * alloc, Node ** root);
static int yylex (Node ** lvalp, char ** input, Allocator * alloc);
static void yyerror (char ** input, Allocator * alloc,
		     Node ** root, char const *);
%}

%param {char ** input}
%param {Allocator * alloc}
%parse-param {Node ** root}
%define api.pure full
%define api.value.type {Node *}

%token NUM        /* Simple double precision number.  */
%token VAR FNCT   /* Variable and Function.  */

%left '?' ':'
%left OR
%left AND
%left EQUAL DIFF
%left '<' LE '>' GE
%left '-' '+'
%left '*' '/'
%left NEG '!'    /* negation--unary minus */
%right '^'    /* exponentiation */

%%
input:
exp      { *root = $1;   }
| error  { *root = NULL; YYABORT; }
;

exp:
NUM

| VAR

| VAR '[' ']'

| VAR '[' exp ']'    {
  $$ = $1; $1->e[0] = $3;
}

| VAR '[' exp ',' exp ']' {
  $$ = $1; $1->e[0] = $3; $1->e[1] = $5;
}

| VAR '[' exp ',' exp ',' exp ']' {
  $$ = $1; $1->e[0] = $3; $1->e[1] = $5; $1->e[2] = $7;
}

| FNCT '(' exp ')'   {
  $$ = $1; $1->e[0] = $3;
}

| exp '+' exp        {
  $$ = new_expr (alloc, '+', $1, $3, NULL);
}

| exp '-' exp        {
  $$ = new_expr (alloc, '-', $1, $3, NULL);
}

| exp '*' exp        {
  $$ = new_expr (alloc, '*', $1, $3, NULL);
}

| exp '/' exp        {
  $$ = new_expr (alloc, '/', $1, $3, NULL);
}

| '-' exp  %prec NEG {
  $$ = new_expr (alloc, 'm', $2, NULL, NULL);
}

| '+' exp  %prec NEG

| exp '^' exp        {
  $$ = new_expr (alloc, '^', $1, $3, NULL);
}

| exp '>' exp        {
  $$ = new_expr (alloc, '>', $1, $3, NULL);
}

| exp '<' exp        {
  $$ = new_expr (alloc, '<', $1, $3, NULL);
}

| exp LE exp        {
  $$ = new_expr (alloc, 'L', $1, $3, NULL);
}

| exp GE exp        {
  $$ = new_expr (alloc, 'G', $1, $3, NULL);
}

| exp EQUAL exp        {
  $$ = new_expr (alloc, '=', $1, $3, NULL);
}

| exp DIFF exp        {
  $$ = new_expr (alloc, 'i', $1, $3, NULL);
}

| exp OR exp        {
  $$ = new_expr (alloc, 'O', $1, $3, NULL);
}

| exp AND exp        {
  $$ = new_expr (alloc, 'A', $1, $3, NULL);
}

| exp '?' exp ':' exp {
  $$ = new_expr (alloc, '?', $1, $3, $5);
}

| '(' exp ')'        {
  $$ = $2;
}
;
%%

/* Called by yyparse on error.  */
static void
yyerror (char ** input, Allocator * alloc, Node ** root, char const *s)
{
  // printf ("%s\n", s);
}

struct init
{
  char const *fname;
  double (*fnct) (double);
};

static struct init arith_fncts[] =
{
  {"abs",   fabs},
  {"fabs",  fabs},
  {"exp",   exp},
  {"log",   log},
  {"ln",    log},
  {"log10", log10},
  {"sqrt",  sqrt},
  {"sin",   sin},
  {"cos",   cos},
  {"tan",   tan},
  {"asin",  asin},
  {"acos",  acos},
  {"atan",  atan},
  {"sinh",  sinh},
  {"cosh",  cosh},
  {"tanh",  tanh},
  {"asinh", asinh},
  {"acosh", acosh},
  {"atanh", atanh},
  {"erf",   erf},
  {"erfc",  erfc},
  {NULL, NULL}
};

static int
yylex (YYSTYPE *lvalp, char ** input, Allocator * alloc)
{
  int c;

  /* Ignore white space, get first nonwhite character.  */
  while ((c = *((*input)++)) == ' ' || c == '\t' || c == '\n' || c == '\t');

  if (c == '\0')
    return EOF;

  /* Char starts a number => parse the number.         */
  if (c == '.' || isdigit (c))
    {
      (*input)--;
      double val;
      int pos;
      if (sscanf (*input, "%lf%n", &val, &pos) == 1) {
	*lvalp = new_value (alloc, val);
	*input += pos;
      }
      else
	return EOF;
      return NUM;
    }

  /* Char starts an identifier => read the name.       */
  if (isalpha (c))
    {
      /* Initially make the buffer long enough
         for a 40-character symbol name.  */
      int length = 40;
      char * symbuf = (char *)malloc (length + 1);
      int i = 0;
      do
        {
          /* If buffer is full, make it bigger.        */
          if (i == length)
            {
              length *= 2;
              symbuf = (char *) realloc (symbuf, length + 1);
            }
          /* Add this character to the buffer.         */
          symbuf[i++] = c;
          /* Get another character.                    */
          c = *((*input)++);
        }
      while (isalnum (c) || c == '_' || c == '.');

      (*input)--;
      symbuf[i] = '\0';

      Node * n = allocate (alloc, sizeof (Node));
      
      struct init * p = arith_fncts;      
      int type = VAR;
      while (p->fname) {
	if (!strcmp (symbuf, p->fname)) {
	  n->type = 'f';
	  n->d.func = p->fnct;
	  type = FNCT;
	  break;
	}
	p++;
      }
      if (type == VAR) {
	n->type = 'v';
	n->d.id = allocate (alloc, strlen(symbuf) + 1);
	memcpy (n->d.id, symbuf, strlen(symbuf) + 1);
      }
      free (symbuf);
      *lvalp = n;
      return type;
    }

  if (c == '|' && **input == '|') {
    (*input)++;
    return OR;
  }

  if (c == '&' && **input == '&') {
    (*input)++;
    return AND;
  }

  if (c == '=' && **input == '=') {
    (*input)++;
    return EQUAL;
  }  

  if (c == '!' && **input == '=') {
    (*input)++;
    return DIFF;
  }

  if (c == '>' && **input == '=') {
    (*input)++;
    return GE;
  }

  if (c == '<' && **input == '=') {
    (*input)++;
    return LE;
  }
  
  /* Any other character is a token by itself.        */
  return c;
}

Node * copy_node (Node * n)
{
  Node * a = malloc (sizeof (Node));
  memcpy (a, n, sizeof (Node));
  if (a->type == 'v')
    a->d.id = strdup (n->d.id);
  int i;
  for (i = 0; i < 3; i++)
    if (a->e[i])
      a->e[i] = copy_node (a->e[i]);
  return a;
}

Node * parse_node (char * code)
{
  Node * root = NULL;
  Allocator * alloc = new_allocator();
  int status = yyparse (&code, alloc, &root);
  if (status == 0)
    root = copy_node (root);
  else
    root = NULL;
  free_allocator (alloc);
  return root;
}

void free_node (Node * n)
{
  int i;
  for (i = 0; i < 3; i++)
    if (n->e[i])
      free_node (n->e[i]);
  if (n->type == 'v')
    free (n->d.id);
  free (n);
}

void print_node (Node * n, FILE * fp)
{
  fprintf (fp, "n%p ", n);
  if (n->type == '1')
    fprintf (fp, "[label=\"%g\"];\n", n->d.value);
  else if (n->type == 'v')
    fprintf (fp, "[label=\"%s\"];\n", n->d.id);
  else if (n->type == 'f')
    fprintf (fp, "[label=\"%p\"];\n", n->d.func);
  else
    fprintf (fp, "[label=\"%c\", shape=box];\n", n->type);
  int i;
  for (i = 0; i < 3; i++)
    if (n->e[i])
      fprintf (fp, "n%p -> n%p;\n", n, n->e[i]);
  for (i = 0; i < 3; i++)
    if (n->e[i])
      print_node (n->e[i], fp);
}

void reset_node_type (Node * n, char type)
{
  n->type = type; free (n->d.id);
}

#if STANDALONE
int
main (int argc, char * argv[])
{
  assert (argc > 1);
  Node * root = parse_node (argv[1]);
  if (root) {
    printf ("digraph G {\n");
    print_node (root, stdout);
    printf ("}\n");
    free_node (root);
  }
  else
    fprintf (stderr, "syntax error\n");
}
#endif // STANDALONE
