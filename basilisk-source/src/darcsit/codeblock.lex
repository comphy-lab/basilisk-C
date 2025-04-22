%option reentrant noyywrap extra-type="struct MyScanner *"
%{
  #include <stdlib.h>
   
  #undef YY_BUF_SIZE
  #define YY_BUF_SIZE 262144

  typedef struct {
    char id[80], file[512], line[80];
  } Tag;

  struct MyScanner {
    char * page;
    Tag * decl, * call, * incl;
    int ndecl, ncall, nincl;
    int i, intypedef, intypedefscope, scope;
    int incode;
  };

  typedef struct { char * s, * sub; } Subst;

  #define lookup_decl(x) lookup_tag(x, yyextra->decl, yyextra->ndecl);
  #define lookup_call(x) lookup_tag(x, yyextra->call, yyextra->ncall);
  static Tag * lookup_tag (const char * id, Tag * tags, int ntags) {
    int i;
    for (i = 0; i < ntags; i++)
      if (!strcmp (id, tags[i].id))
	return &tags[i];
    return NULL;
  }

  #define output_c(c) putchar(c)
  #define echo() output_s(yytext)
  #define output_s(s) fputs(s, stdout)

  static char * baseurl = "", * ext = NULL, * basilisk_url = NULL;

  char * url (char * s, char * baseurl)
  {
    static char s1[256];
    if (s[0] == '/') {
      // assert (strlen(baseurl) + strlen(s) + 1) < 256;
      char * root = getenv ("DOCUMENT_ROOT");
      if (root && s[strlen(root)] == '/' && !strncmp (root, s, strlen(root)))
	s += strlen(root);
      strcat (strcpy (s1, baseurl), s);
    }
    else
      strcpy (s1, s);
    if (ext) {
      char page[strlen(s1) + strlen (".page") + 1];
      strcpy (page, s1);
      char * anchor = strchr (page, '#');
      if (anchor)
	*anchor = '\0';
      strcat (page, ".page");
      // fprintf (stderr, "testing '%s'\n", page);
      FILE * fp = fopen (page, "r");
      if (fp) {
	if (anchor) {
	  strcpy (page, s1);
	  *strchr (s1, '#') = '\0';
	  strcat (s1, ".html");
	  strcat (s1, anchor);
	}
	else
	  strcat (s1, ".html");
	fclose (fp);
      }
    }
    return s1;
  }

# define INCODE() (yyextra->incode)
  
  static int check_tag (char * text, struct MyScanner * scan) {
    Tag * t = lookup_tag (text, scan->call, scan->ncall);
    if (t) {
      // link to another page
      output_s ("<a href=");
      if (!strcmp(t->file, "stdlib")) {
	output_s ("http://man7.org/linux/man-pages/man3/");
	output_s (text);
	output_s (".3.html");
      }
      else if (!strcmp(t->file, "basilisk")) {
	output_s (url ("/Basilisk%20C", basilisk_url));
	output_c ('#');
	output_s (t->line);
      }
      else {
	output_s (url (t->file, basilisk_url));
	output_c ('#');
	char s1[256];
	if (t->file[0] == '/')
	  strcat (strcpy (s1, baseurl), t->file);
	else
	  strcpy (s1, t->file);
	output_s (t->line);
      }
      output_s (">");
      output_s (text);
      output_s ("</a" ">");
    }
    return (t != NULL);
  }

  static char * append_c (char * tmp, int c) {
    int len = strlen(tmp);
    tmp = realloc (tmp, len + 2);
    tmp[len] = c; tmp[len+1] = '\0';
    return tmp;
  }

  static char * append_s (char * tmp, char * s) {
    while (*s)
      tmp = append_c (tmp, *s++);
    return tmp;
  }

  #define nonspace(s) { while (strchr(" \t\v\n\f", *s)) s++; }
  #define space(s) { while (!strchr(" \t\v\n\f", *s)) s++; }
  static void comment(yyscan_t scanner);
%}

SYMID  [a-zA-Z0-9]
ID  [a-zA-Z0-9_]
SP  [ \t]
ES  (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
WS  [ \t\v\n\f]

%%

\v[^\v]*\v {
  // code line numbers
  output_s ("<span id=");
  char * s = yytext;
  while (*s) {
    if (strchr ("0123456789", *s))
      output_c (*s);
    s++;
  }
  output_s ("></span>");
}

^{SP}*"@def"{SP} {
  yyextra->incode = 0; echo();
}

^{SP}*"@"{SP}* {
  yyextra->incode = 1; echo();
}

"<code class=\"sourceCode c\">"  {
  yyextra->incode = 1; echo();
}

"</code>"                        {
  yyextra->incode = 0; echo();
}

\&quot;({ID}|[-/])+\.h\&quot; |
">"#{SP}*include{SP}+\".*\""<" |
^{SP}*#{SP}*include{SP}+\".*\" {
  if (!INCODE())
    REJECT;
  // add HTML links in include headers
  char c = '"', * s = strchr(yytext, c);
  if (!s)
    c = ';', s = strchr(yytext, c);
  *s++ = '\0';
  char c1 = '"', * s1 = strchr(s, c1);
  if (!s1)
    c1 = '&', s1 = strchr(s, c1);
  *s1++ = '\0';
  // look for header in tags
  char * header = NULL;
  int i = 0;
  for (i = 0; i < yyextra->nincl && !header; i++)
    if (strstr(yyextra->incl[i].id, s))
      header = yyextra->incl[i].id;
  echo();
  output_c (c);
  if (header) {
    output_s ("<a href=");
    output_s (url (header, basilisk_url));
    output_s (">");
  }
  output_s (s);
  if (header)
    output_s ("</a" ">");
  else
    fprintf (stderr, 
	     "codeblock: warning: %s: tag for \"%s\" not found\n",
	     yyextra->page, s);
  output_c (c1);
  output_s (s1);
}

{ID}+{SP}*\( {
  if (!INCODE() || yyextra->scope > 0)
    REJECT;
  // keyword  anchor (function definition)
  // this regexp needs to match that in src/include.lex (for tags)
  char * s = strdup(yytext), * s1 = s;
  while (!strchr(" \t\v\n\f(", *s1)) s1++;
  int c = *s1;
  *s1++ = '\0';
  char * tmp = malloc (sizeof(char)); *tmp = '\0';
  int p = 0, para = 1, c1;
  while (para > p && (c1 = input(yyscanner)) > 0) {
    if (c1 == '\v') {
      // line number
      tmp = append_s (tmp, "<span id=");
      while ((c1 = input(yyscanner)) > 0 && c1 != '\v')
	if (strchr ("0123456789", c1))
	  tmp = append_c (tmp, c1);
      tmp = append_s (tmp, "></span>");
    }
    else {
      tmp = append_c (tmp, c1);
      if (c1 == '(') para++;
      else if (c1 == ')') para--;
    }
  }
  if (c1 != ')') {
    output_s (s);
  }
  else {
    while ((c1 = input(yyscanner)) > 0) {
      if (c1 == '\v') {
	// line number
	tmp = append_s (tmp, "<span id=");
	while ((c1 = input(yyscanner)) > 0 && c1 != '\v')
	  if (strchr ("0123456789", c1))
	    tmp = append_c (tmp, c1);
	tmp = append_s (tmp, "></span>");
      }
      else {
	tmp = append_c (tmp, c1);
	if (c1 == '{' || c1 == ';')
	  break;
	if (!strchr(" \t\v\n\f", c1))
	  break;
      }
    }
    if (c1 != '{') {
      output_s (s);
    }
    else {
      yyextra->scope++;
      Tag * t = lookup_decl (s);
      if (t == NULL || !strstr (yyextra->page, t->file)) {
	if (t != NULL)
	  fprintf (stderr, "file: %s page: %s\n", t->file, yyextra->page);
#if 1
	fprintf (stderr, 
		 "codeblock: warning: %s: tag for '%s' not found\n",
		 yyextra->page, s);
#endif
	output_s (s);
      }
      else
	printf ("<a id=%s>%s</a>", t->id, s);
    }
  }
  output_c (c);
  output_s (s1);
  output_s (tmp);
  free (tmp);
  free (s);
}

\'.\' {
  echo(); // quoted character
}

"\&#39;"."\&#39;" {
  echo(); // quoted character in HTML
}

"{" {
  echo();
  yyextra->scope++;
}

"}" {
  echo();
  yyextra->scope--;
  if (yyextra->scope < 0) {
    fprintf (stderr, "warning: %s: error: mismatched '}'\n", yyextra->page);
    yyextra->scope = 0;
  }
}

"<img src=\""[^\"]+\" |
"<a href=\""[^\"]+\" {
  // URL
  yytext[strlen(yytext) - 1] = '\0';
  char * u = strchr (yytext, '"');
  *u = '\0';
  output_s (yytext);
  output_c ('"');
  output_s (url (u + 1, baseurl));
  output_c ('"');
}

typedef {
  if (!INCODE())
    REJECT;
  echo();
  yyextra->intypedef = 1; yyextra->intypedefscope = yyextra->scope;
}

{ID}+{WS}*; {
  if (!INCODE())
    REJECT;
  // keyword  anchor (typedef)
  // this regexp needs to match that in src/include.lex (for tags)
  if (yyextra->intypedef && yyextra->intypedefscope == yyextra->scope) {
    char * s = yytext; space(s); *s-- = '\0';
    if (*s == ';')
      *s = '\0';
    yyextra->intypedef = 0;
    s = yytext;
    Tag * t = lookup_decl (s);
    if (t == NULL || !strstr (yyextra->page, t->file)) {
      if (t != NULL)
	fprintf (stderr, "file: %s page: %s\n", t->file, yyextra->page);
      fprintf (stderr, 
	       "codeblock: warning: %s: tag for '%s' not found\n",
	       yyextra->page, s);
      output_s (s);
    }
    else
      printf ("<a id=%s>%s</a>", t->id, s);
    output_c (';');
  }
  else
    REJECT;  
}

scalar|vector|tensor {
  output_s ("<span class=\"dt\">");
  echo();
  output_s ("</span>");
}

{ID}+ {
  if (!INCODE())
    REJECT;
  // keyword links
  if (!check_tag (yytext, yyextra))
    REJECT;
}

"/*" {
  if (!INCODE())
    REJECT;
  echo();
  comment(yyscanner);
}

{ID}+                                 echo();
.                                     echo();
\n                                    echo();
\"([^\"\\\n]|{ES})*\"  { echo(); /* STRING_LITERAL */ }

%%

static void comment(yyscan_t scanner)
{
  int c;
  while ((c = input(scanner)) > 0) {
    if (c == '*') {
      output_c (c);
      while ((c = input(scanner)) == '*')
	output_c (c);
      output_c (c);
      if (c == '/')
	return;
      if (c == 0)
	break;
    }
    else if (c == '\v') {
      // line number
      output_s ("<span id=");
      while ((c = input(scanner)) && strchr ("0123456789", c))
	output_c (c);
      output_s ("></span>");
    }
    else
      output_c (c);
  }
  fprintf (stderr, "codeblock: warning: %s: unterminated comment\n", 
	   yyget_extra(scanner)->page);
}

static void read_tagfile (struct MyScanner * scan)
{
  char files[2][200] = {"", ""};
  char * s = getenv ("BASILISK");
  if (s)
    strcpy (files[0], s);
  strcat (files[0], "/external.tags");
  strcpy (files[1], scan->page);
  strcat (files[1], ".tags");
  int i;
  for (i = 0; i < 2; i++) {
    char * s = files[i];
    FILE * fp = fopen (s, "r");
    if (fp == NULL)
      perror (s);
    else {
      Tag t;
      char type[10];
      while (fscanf (fp, "%s %s %s %s", type, t.id, t.file, t.line) == 4) {
	if (!strcmp (type, "decl")) {
	  scan->ndecl++;
	  scan->decl = realloc (scan->decl, scan->ndecl*sizeof(Tag));
	  scan->decl[scan->ndecl-1] = t;
	}
	else if (!strcmp (type, "call")) {
	  scan->ncall++;
	  scan->call = realloc (scan->call, scan->ncall*sizeof(Tag));
	  scan->call[scan->ncall-1] = t;
	}
	else if (!strcmp (type, "incl")) {
	  scan->nincl++;
	  scan->incl = realloc (scan->incl, scan->nincl*sizeof(Tag));
	  scan->incl[scan->nincl-1] = t;
	}
      }
      fclose (fp);
    }
  }
}

void codeblock (char * page)
{
  yyscan_t scanner;
  struct MyScanner sdata = {};
  sdata.i = sdata.scope = sdata.intypedef = 0;
  sdata.page = page;
  sdata.incode = 0;
  read_tagfile (&sdata);

  yylex_init_extra (&sdata, &scanner);
  yyset_in (stdin, scanner);
  yyset_out (stdout, scanner);
  yyset_debug (2, scanner);
  yylex (scanner);
  yylex_destroy (scanner);
  free (sdata.decl);
  free (sdata.call);
  free (sdata.incl);
}

static void usage () {
  fprintf (stderr,
	   "usage: codeblock BASEURL FILE.[ch] [EXT] < FILE.[ch].html\n");
  exit (1);
}

int main (int argc, char * argv[])
{
  if (argc < 3)
    usage();
  char * name = argv[2];

  baseurl = argv[1];
  if (argc >= 4)
    ext = argv[3];
  basilisk_url = getenv ("HTTP_BASILISK_URL");
  if (basilisk_url == NULL)
    basilisk_url = baseurl;
  
  codeblock (name);
  return 0;
}
