%option noyywrap
%option yylineno
%{
  #include <assert.h>
  #include <sys/stat.h>
  #include <sys/types.h>
  #include "ast/allocator.h"
  
  enum { FUNCTION, TYPEDEF };

  typedef struct {
    char * id, * file;
    int line, type;
  } Tag;

  static Tag * tagsa = NULL;
  static int ntags = 0, target = 1, keywords_only = 0, scope = 0, intypedef = 0;
  static int warninclude = 0;
  static FILE * swigfp = NULL;
  static Allocator * alloc = NULL;
  
  static void append_tag (Tag t) {
    ntags++;
    tagsa = realloc (tagsa, ntags*sizeof(Tag));
    tagsa[ntags-1] = t;
    tagsa[ntags-1].id = strdup (t.id);
    tagsa[ntags-1].file = strdup (t.file);
    tagsa[ntags-1].type = t.type;
    char * page = strstr (tagsa[ntags-1].file, ".page");
    if (page) *page = '\0';
  }

  static Tag * lookup_tag (const char * id) {
    int i;
    for (i = 0; i < ntags; i++)
      if (!strcmp(tagsa[i].id, id))
	return &tagsa[i];
    return NULL;
  }

  static FILE * fdepend = NULL, * ftags = NULL, * myout = NULL;
  static char * fname;
  
  static char * paths[100] = { LIBDIR }, grid[80] = "";
  static int npath = 1, hasgrid = 0, debug = 0;
  static int dimension = 0, bghosts = 0, layers = 0, gpu = 0;
  static int incode;    // are we in code (or in a code block)?
  
  static char * strip_path (char * s) {
    char * s1 = s;
    do {
      while (*s1 != '/' && *s1 != '\0') s1++;
      if (*s1 == '\0') return s;
    } while (*(++s1) != '/');
    while (*s1 == '/') s1++;
    return s1;
  }
  
  static char * _stack[100]; int stack = -1;

  static char * push (char * s) {
    assert (stack + 1 < 100);
    char * f = allocate (alloc, strlen (s) + 1);
    strcpy (f, s);
    _stack[++stack] = f;
    return f;    
  }

  static int strcmps (char * s1, char * s2) {
    while (*s1 != '\0' && *s2 != '\0') {
      if (*s1 != *s2)
	return 1;
      if (*s1 == '/') {
	while (*s1 == '/') s1++;
	while (*s2 == '/') s2++;
      }
      else
	s1++, s2++;
    }
    return !(*s1 == '\0' && *s2 == '\0');
  }

  static void push_once (char * s) {
    static char * _processed[100];
    static int processed = 0;
    int i;
    char * s1 = strip_path (s);
    for (i = 0; i < processed; i++)
      if (!strcmps (s1, strip_path (_processed[i])))
	return; // already processed    
    assert (processed < 100);
    _processed[processed++] = push (s);
  }
  
  #define pop()  _stack[stack--];

  static void singleslash (char * path, FILE * fp)
  {
    char * s = path, slash = 0;
    while (*s != '\0') {
      if (*s == '/') {
	if (!slash)
	  fputc (*s, fp);
	slash = 1;
      }
      else {
	slash = 0;
	fputc (*s, fp);
      }
      s++;
    }
  }

  static FILE * openpath (const char * name, const char * mode, char ** path) {
    int i;
    for (i = npath; i >= 0; i--) {
      char * p = allocate (alloc, strlen (paths[i]) + strlen (name) + 3);
      strcpy (p, paths[i]); strcat (p, "//"); strcat (p, name);
      FILE * fp = fopen (p, mode);
      if (fp) {
	if (fdepend) {
	  fputc ('\t', fdepend); singleslash (p, fdepend);
	  fputs (" \\\n", fdepend);
	}
	*path = p;
	return fp;
      }
    }
    return NULL;
  }

#define nonspace(s) { while (strchr(" \t\v\n\f", *s)) s++; }
#define space(s) { while (!strchr(" \t\v\n\f", *s)) s++; }

  static char * shortpath (char * path) {
    char * file = strstr (path, LIBDIR);
    if (file == path)
      return file + strlen(LIBDIR) - strlen("src") - 1; // remove root
    else
      return path;
  }

  static void check_tag (char * text) {
    if (ftags && keywords_only && incode) {
      Tag * t;
      if (target && (t = lookup_tag(text))) {
	switch (t->type) {
	case FUNCTION:
	  if (debug)
	    fprintf (stderr, "%s:%d: function call '%s'\n", 
		     fname, yylineno, text);
	  break;
	case TYPEDEF:
	  if (debug)
	    fprintf (stderr, "%s:%d: typedef reference '%s'\n", 
		     fname, yylineno, text);
	  break;
	}
	fprintf (ftags, "call %s %s %s\n", 
		 t->id, shortpath (t->file), t->id);
      }
    }
  }

  static int yyerror(const char * s);
  static int comment(void);
  static void echo() {
    if (myout) {
      if (incode)
	fputs (yytext, myout);
      else { // only keep newlines
	char * s = yytext;
	while (*s != '\0') {
	  if (*s == '\n')
	    fputc ('\n', myout);
	  s++;
	}
      }
    }
  }

  static void echo_c (int c) {
    if (myout) {
      if (incode)
	fputc (c, myout);
      else if (c == '\n') // only keep newlines
	fputc ('\n', myout);
    }
  }
%}

ID     [a-zA-Z0-9_]
SP     [ \t]
WS     [ \t\v\n\f]
ES     (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
BEGINCODE ^[SP]*[~]{3,}(c|literatec)[^\n]*\n
ENDCODE   ^[SP]*[~]{3,}[^\n]*\n
FDECL     {ID}+{SP}*\(

%%

{BEGINCODE} {
  if (incode) {
    yylineno--;
    return yyerror ("code blocks cannot be nested");
  }
  incode = 1;
  if (myout) fputc ('\n', myout);
}

{ENDCODE} {
  if (incode) {
    incode = 0;
    if (myout) fputc ('\n', myout);
  }
}

\'.\' {
  echo(); // quoted character
}

"{" {
  scope++;
}

"}" {
  scope--;
  if (scope < 0) {
    if (warninclude)
      fprintf (stderr, "%s:%d: warning: mismatched '}'\n", fname, yylineno);
    else {
      fprintf (stderr, "%s:%d: error: mismatched '}'\n", fname, yylineno);
      exit (1);
    }
    scope = 0;
  }
}

^{SP}*#{SP}*include{SP}+\"[^\"]*\"[^\n]*\n {
  // include "..."
  if (fdepend && strstr (yytext, "// nodep"))
    return 0;
  echo();
  if (!keywords_only) {
    char * s = strchr(yytext, '"');
    s++;
    char * e = &s[strlen(s) - 1];
    while (*e != '"') {
      *e = '\0'; e--;
    }
    *e = '\0';
    char * path;
    FILE * fp = openpath (s, "r", &path);
    if (fp != NULL) {
      push_once (path);
      if (ftags && target) {
	fputs ("incl ", ftags);
	singleslash (shortpath(path), ftags);
	fprintf (ftags, " %s %d\n", fname, yylineno - 1);
      }
      if (swigfp && target) {
	char * dot = strstr (path, ".h");
	if (dot) {
	  strcpy (dot, ".i");
	  fputs ("%include \"", swigfp);
	  singleslash (path, swigfp);
	  fputs ("\"\n", swigfp);
	}
      }
      fclose (fp);
    }
    else {
      fprintf (stderr, "%s:%d: %s: %s: No such file or directory\n", 
	       fname, yylineno - 1, warninclude ? "warning" : "error", s);
      if (!warninclude)
	return 1;
    }
  }
}

^{SP}*#{SP}*define{SP}+GRIDNAME{WS}+ {
    echo();
    if (grid[0] == '\0' && !hasgrid) {
      hasgrid = 1;
      char * s = fname;
      while (strchr (s, '/')) {
	s = strchr (s, '/') + 1;
	if (!strncmp (s, "grid/", 5)) {
	  s += 5;
	  break;
	}
      }
      strcpy (grid, s);
      if ((s = strchr (grid, '.'))) *s = '\0';
    }
}

^{SP}*#{SP}*define{SP}+dimension{WS}+[123]{SP}*$ {
  char * s = strstr (yytext, "dimension");
  space(s); nonspace(s);
  dimension = atoi(s);
}

^{SP}*#{SP}*define{SP}+BGHOSTS{WS}+[12]{SP}*$ {
  char * s = strstr (yytext, "BGHOSTS");
  space(s); nonspace(s);
  bghosts = atoi(s);
}

^{SP}*#{SP}*define{SP}+_GPU{WS}+1{SP}*$ {
  gpu = 1;
}

^{SP}*#{SP}*define{SP}+LAYERS{WS}+1{SP}*$ {
  layers = 1;
}

^{SP}*{ID}+{SP}*\**({SP}+{ID}+{SP}*\**)*{SP}+{ID}+{SP}*\( {
  // function definition
  echo();
  if (ftags && scope == 0) {
    //    fprintf (stderr, "'%s'\n", yytext);
    char * s = yytext; int nl = 0;
    int fstatic = 0;
    s = strtok (s, " \t\v\n\f(");
    char * id = s;
    while (s) {
#if 0 // ignore static functions
      if (!strcmp (s, "static"))
	fstatic = 1;
#endif
      id = s;
      s = strtok (NULL, " \t\v\n\f(");
      if (s)
	check_tag (id);
    }
    s = strdup (id);
    if (!fstatic && !keywords_only && strcmp(s, "if")) {
      //      fprintf (stderr, "id: '%s'\n", s);
      Tag t = { s, fname, yylineno - nl, FUNCTION};
      int p = 0, para = 1, c;
      while (para > p && (c = input())) {
	echo_c (c);
	if (c == '(') para++;
	else if (c == ')') para--;
      }
      if (c == ')') {
	while ((c = input())) {
	  echo_c (c);
	  if (c == '{' || c == ';')
	    break;
	  if (!strchr(" \t\v\n\f", c))
	    break;
	}
	if (c == '{') {
	  scope++;
	  append_tag (t);
	  if (debug)
	    fprintf (stderr, "%s:%d: function declaration '%s'\n", 
		     tagsa[ntags-1].file, tagsa[ntags-1].line, 
		     tagsa[ntags-1].id);
	  if (target)
	    fprintf (ftags, "decl %s %s %d\n", 
		     tagsa[ntags-1].id, tagsa[ntags-1].file, 
		     tagsa[ntags-1].line);  
	}
      }
    }
    free (s);
  }
}

typedef{WS}+ {
  echo();
  if (ftags && !keywords_only)
    intypedef = scope + 1;
}

{ID}+{WS}*; {
  if (intypedef && scope == intypedef - 1) {
    echo();
    char * s = yytext; space(s); *s-- = '\0';
    if (*s == ';')
      *s = '\0';
    Tag t = { yytext, fname, yylineno, TYPEDEF};
    append_tag (t);
    if (debug)
      fprintf (stderr, "%s:%d: typedef '%s'\n", 
	       tagsa[ntags-1].file, tagsa[ntags-1].line, tagsa[ntags-1].id);
    if (target)
      fprintf (ftags, "decl %s %s %d\n", 
	       tagsa[ntags-1].id, tagsa[ntags-1].file, tagsa[ntags-1].line);
    intypedef = 0;
  }
  else
    REJECT;  
}

{ID}+ {
  // keyword in target
  echo();
  check_tag (yytext);
}

"/*"              { echo(); if (incode && comment()) return 1; }
"//".*            {
  if (!incode)
    REJECT;
  /* consume //-comment */
  echo();
}
.                   echo();
[\n]                echo();
({SP}?\"([^\"\\\n]|{ES})*\"{WS}*)+  echo(); /* STRING_LITERAL */

%%

int yyerror (const char * s)
{
  if (0) yyunput (0, NULL); // just prevents 'yyunput unused' compiler warning
  fprintf (stderr, "%s:%d: error: %s\n", fname, yylineno, s);
  return 1;
}

static int getput(void)
{
  int c = input();
  if (myout)
    fputc (c, myout);
  return c;
}

static int comment(void)
{
  int c, lineno = yylineno;
  while ((c = getput())) {
    if (c == '*') {
      while ((c = getput()) == '*')
	;
      if (c == '/')
	return 0;
      if (c == 0)
	break;
    }
  }
  fprintf (stderr, "%s:%d: warning: unterminated comment\n", fname, lineno);
  return 0;
}

void stripname (char * path)
{
  char * s = &path[strlen(path)];
  while (s != path && *s != '/')
    *s-- = '\0';
  if (s == path)
    strcpy (path, ".");
  else
    *s = '\0';
}

char * stripslash (char * path)
{
  char * strip = malloc (strlen (path) + 1), * s = path, * o = strip;
  int slash = 0;
  do {
    if (*s == '/') {
      if (!slash)
	*o++ = *s;
      slash = 1;
    }
    else {
      *o++ = *s;
      slash = 0;
    }
  } while (*s++ != '\0');
  return strip;
}

static int is_code (const char * file)
{
  // check whether file has a .c or .h extension
  char * s = strstr (file, ".c");
  if (!s)
    s = strstr (file, ".h");
  return s && (s[2] == '\0' || s[2] == '.');
}

static int include (char * file, FILE * fin, FILE * fout)
{
  fname = stripslash (file);
  paths[npath] = allocate (alloc, strlen (file) + 1);
  strcpy (paths[npath], file);
  stripname (paths[npath]);
  yyin = fin;
  myout = fout;
  yylineno = 1;
  scope = intypedef = 0;
  long header = fout ? ftell (fout) : 0;
  incode = is_code (file);
  //  yy_flex_debug = 1;
  int ret = yylex();
  yylex_destroy();
  if (fout && incode) {
    // Assume the entire file is pure code
    fseek (fout, header, SEEK_SET);
    rewind (fin);
    char s[81];
    while (fgets (s, 81, fin))
      fputs (s, fout);
  }
  free (fname);
  return ret;
}

FILE * writepath (char * path, const char * mode);

static void compdir (char * file, const char * dir)
{
  push (file);
  while (stack >= 0) {
    char * path = pop();
    FILE * fin = fopen (path, "r");
    if (fin == NULL) {
      perror (path);
      exit (1);
    }
    FILE * fout = NULL;
    if (dir) {
      char * file = strstr (path, "//");
      if (file) file += 2; else file = path;
      char * out = malloc (strlen (dir) + strlen (file) + 2);
      strcpy (out, dir);
      strcat (out, "/");
      strcat (out, file);
      fout = writepath (out, "w");
      if (fout == NULL) {
	perror (out);
	exit (1);
      }
      free (out);

      // all headers are included once only
      static int nf = 0;
      fprintf (fout, 
	       "#ifndef BASILISK_HEADER_%d\n"
	       "#define BASILISK_HEADER_%d\n", nf, nf);
      nf++;

      fputs ("#line 1 \"", fout);
      singleslash (path, fout);
      fputs ("\"\n", fout);
    }
    if (include (path, fin, fout))
      exit (1);
    fclose (fin);
    if (fout) {
      fputs ("\n#endif\n", fout);
      fclose (fout);
    }
    target = 0;
  }
}

static void prepend_path (char * path)
{
  int j;
  for (j = npath; j > 0; j--)
    paths[j] = paths[j-1];
  paths[0] = path;
  npath++;
}

void includes (int argc, char ** argv,
	       char ** grid1, int * default_grid,
	       int * dim, int * bg, int * lyrs, int * gpus,
	       const char * dir)
{
  int depend = 0, tags = 0, swig = 0;
  char * file = NULL, * output = NULL;
  int i;
  warninclude = 0;
  alloc = new_allocator();
  char * basilisk_include_path = getenv ("BASILISK_INCLUDE_PATH");
  if (basilisk_include_path) {
    basilisk_include_path = strdup (basilisk_include_path);
    char * s = strtok (basilisk_include_path, ":");
    while (s) {
      prepend_path (s);
      s = strtok (NULL, ":");
    }
  }
  for (i = 1; i < argc; i++) {
    if (!strncmp (argv[i], "-grid=", 6))
      strcpy (grid, &argv[i][6]);
    else if (!strcmp (argv[i], "-MD"))
      depend = warninclude = 1;
    else if (!strcmp (argv[i], "-tags"))
      tags = warninclude = 1;
    else if (!strcmp (argv[i], "-python"))
      swig = 1;
    else if (!strcmp (argv[i], "-debug"))
      debug = 1;
    else if (!strcmp (argv[i], "-o"))
      output = argv[++i];
    else if (!strncmp (argv[i], "-I", 2))
      prepend_path (argv[i] + 2);
    else if (argv[i][0] != '-' && \
	     (tags || !strcmp (&argv[i][strlen(argv[i]) - 2], ".c"))) {
      if (file) {
	fprintf (stderr, "usage: include [OPTIONS] FILE.c\n");
	exit (1);
      }
      file = argv[i];
    }
  }
  if (depend && file) {
    if (!output) output = file;
    char ndep[80], * s = &output[strlen(output)-1];
    while (*s != '.' && s != output) s--;
    if (output != file || s == output)
      /* generate dep files with suffixes included for -o option */
      strcpy (ndep, output);
    else {
      *s = '\0';
      strcpy (ndep, output);
      *s = '.';
    }
    if (strlen(ndep) < 2 || strcmp (&ndep[strlen(ndep)-2], ".d")) {
      if (tags)
	strcat (ndep, ".tags.d");
      else
	strcat (ndep, ".d");
    }
    else
      output[strlen(ndep)-2] = '\0'; // strip trailing ".d";
    fdepend = fopen (ndep, "w");
    if (!fdepend) {
      perror (ndep);
      exit (1);
    }
    char * page = strstr (output, ".page");
    if (tags && page) {
      *page = '\0';
      fprintf (fdepend, "%s.tags:\t\\\n", output);
      *page = '.';
    }
    else
      fprintf (fdepend, "%s:\t\\\n", output);
  }
  else if (tags && file) {
    if (!output) output = file;
    char ndep[80];
    // strip trailing .page
    strcpy (ndep, output);
    char * page = strstr (ndep, ".page");
    if (page)
      *page = '\0';
    strcat (ndep, ".tags");
    ftags = fopen (ndep, "w");
    if (!ftags) {
      perror (ndep);
      exit (1);
    }
  }
  if (file) {
    if (swig) {
      char swigname[80];
      strcpy (swigname, file);
      char * dot = strchr (swigname, '.');
      *dot = '\0'; strcat (swigname, ".i");
      swigfp = fopen (swigname, "w");
      if (!swigfp) {
	fprintf (stderr, "include: could not open '%s': ", swigname);
	perror ("");
	exit (1);
      }
      *dot = '\0';
      fprintf (swigfp, "%%module %s\n", swigname);
      fputs ("%include \"", swigfp);
      fputs (LIBDIR, swigfp);
      fputs ("/common.i\"\n", swigfp);
    }
    target = 1;
    compdir (file, dir);
    if (!hasgrid && is_code (file)) {
      char * path, gridpath[80] = "grid/";
      strcat (gridpath, grid[0] != '\0' ? grid : "quadtree"); strcat (gridpath, ".h");
      FILE * fp = openpath (gridpath, "r", &path);
      if (!fp) {
	fprintf (stderr, "include: invalid grid '%s': ", grid);
	perror ("");
	exit (1);
      }
      fclose (fp);
      target = 0;
      compdir (path, dir);
      hasgrid = 0;
    }
    if (swigfp) {
      char * path, pythonpath[80] = "python.h";
      FILE * fp = openpath (pythonpath, "r", &path);
      if (!fp) {
	perror (pythonpath);
	exit (1);
      }
      fclose (fp);
      target = 0;
      compdir (path, dir);
    }
    if (ftags) {
      // reprocess the target file for keywords
      keywords_only = target = 1;
      compdir (file, dir);
    }
    char * path;    
    FILE * fp = openpath ("common.h", "r", &path);
    assert (fp);
    fclose (fp);
  }
  if (fdepend) {
    fputc ('\n', fdepend);
    fclose (fdepend);
  }
  if (ftags)
    fclose (ftags);
  if (swigfp)
    fclose (swigfp);
  *grid1 = grid;
  *default_grid = !hasgrid;
  if (dimension > 0)
    *dim = dimension;
  *bg = bghosts;
  *lyrs = layers;
  *gpus = gpu;
  free (basilisk_include_path);
  free_allocator (alloc);
}
