%option reentrant noyywrap noinput extra-type="struct MyScanner *"
%{
  #include <stdlib.h>
  #include <stdarg.h>
  #include <sys/wait.h>
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <unistd.h>
  #include <glob.h>
  #include <time.h>
  #include <assert.h>
   
  #undef YY_BUF_SIZE
  #define YY_BUF_SIZE 262144
   
  typedef struct {
    char * error, * warning;
  } Error;
  
  static char C[]      = "\n~~~literatec";
  static char Python[] = "\n~~~python";
  static char Octave[] = "\n~~~matlab";
  static char Bash[] = "\n~~~bash";

  struct MyScanner {
    FILE * in, * out;
    char * page, * basename, * gnuplot, * gnuplot_output, * type, * plotype;
    int i, ncodes, incode, first, line, indent;
    FILE * inbibtex;
    Error * error;
    int nerror, nplots;
  };

#define output_c1(scan, c) fputc(c, scan->out)
#define output_c(c) output_c1(yyextra, c)
#define output_s1(scan, s) fputs(s, scan->out)
#define output_s(s) output_s1(yyextra, s)
  
  static void error_start (struct MyScanner * scan) {
    fprintf (scan->out,
	     "<div class=error id=%d>"
	     "<div id=msg_logo>"
	     "<img src=/img/error.png>"
	     "</div>"
	     "<div id=msg_label>", scan->line + 1);
  }

  static void warning_start (struct MyScanner * scan) {
    fprintf (scan->out,
	     "<div class=message id=%d>"
	     "<div id=msg_logo>"
	     "<img src=/img/warning.png>"
	     "</div>"
	     "<div id=msg_label>", scan->line + 1);
  }

  static void error_end (struct MyScanner * scan) {
    output_s1 (scan, "</div></div>");
  }

  static int check_error (struct MyScanner * scan) {
    int found = 0;
    if (scan->line < scan->nerror) {
      Error * e = &scan->error[scan->line];
      if (e->error || e->warning) {
	if (!scan->first && (scan->incode || scan->plotype))
	  output_s1 (scan, "\n~~~\n");
	if (e->error) {
	  error_start (scan);
	  output_s1 (scan, e->error);
	  free (e->error); e->error = NULL;
	}
	else {
	  warning_start (scan);
	  output_s1 (scan, e->warning);
	  free (e->warning); e->warning = NULL;
	}
	error_end (scan);
	found = 1;
	if (!scan->first && scan->incode && scan->type) {
	  output_s1 (scan, scan->type);
	}
	else if (scan->plotype) {
	  output_c1 (scan, '\n');
	  output_s1 (scan, scan->plotype);
	}
      }
    }
    return found;
  }

  static int spacenb (char * s) {
    int ns = 0;
    while (*s != '\0' && strchr (" \t\v\n\f", *s)) {
      switch (*s) {
      case ' ':  ns++; break;
      case '\t': ns += 8; break;
      case '\n': case '\v': case '\f': ns = 0; break;
      }
      s++;
    }
    return ns;
  }

  static char * acat (const char * s, ...)
  {
    int len = strlen (s) + 1;
    char * c = malloc (len);
    strcpy (c, s);

    va_list ap;
    va_start (ap, s);
    char * s1;
    while ((s1 = va_arg (ap, char *))) {
      len += strlen (s1);
      c = realloc (c, len);
      strcat (c, s1);
    }
    va_end (ap);

    return c;
  }

  static void uline1 (const char * s, struct MyScanner * scan)
  {
    while (*s != '\0')
      if (*s++ == '\n') {
	if (!check_error (scan) && scan->incode && !scan->first &&
	    scan->line > 2)
	  fprintf (scan->out, "\v%d\v", scan->line);
	scan->line++;
      }
  }
  #define uline() uline1(yytext, yyextra)
  
  #define YY_INPUT(buf,result,max_size)			      \
  {							      \
    int c = fgetc (yyextra->in);			      \
    if (c == '\r') /* ignore DOS line endings */	      \
      c = fgetc (yyextra->in);				      \
    if (c == EOF) result = YY_NULL;			      \
    else { buf[0] = c; result = 1; }			      \
  }
%}

ID  [a-zA-Z_0-9]
SP  [ \t]
ES  (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
WS  [ \t\v\n\f]
SCALAR [a-zA-Z_0-9]+[.xyz]*

%%

{WS}*\/[*][*]{SP}* {
  // start of C documentation comment i.e. "/**"
  uline();
  if (yyextra->ncodes > 0 && !yyextra->first)
    output_s ("\n~~~\n");
  output_c ('\n');
  yyextra->incode = 0;
  yyextra->indent = spacenb (yytext);
  yyextra->type = C;
}

{SP}*[*]\/{SP}* {
  uline();
  // end of any C comment block i.e. "*/"
  if (yyextra->incode)
    // end of standard comment
    output_s (yytext);
  else {
    // end of documentation comment
    output_s ("\n");
    yyextra->ncodes++;
    yyextra->incode = 1;
    yyextra->first = 1;
  }
}

{WS}*\"\"\"{SP}* {
  uline();
  if (yyextra->incode) {
    // start of Python documentation comment i.e. """
    if (yyextra->ncodes > 0 && !yyextra->first)
      output_s ("\n~~~\n");
    output_c ('\n');
    yyextra->incode = 0;
    yyextra->indent = spacenb (yytext);
    yyextra->type = Python;
  }
  else {
    // end of Python documentation comment
    output_s ("\n");
    yyextra->ncodes++;
    yyextra->incode = 1;
    yyextra->first = 1;
  }
}

^{WS}*%\{{WS}*$ {
  if (yyextra->incode) {
    uline();
    // start of Octave documentation comment i.e. "%{"
    if (yyextra->ncodes > 0 && !yyextra->first)
      output_s ("\n~~~\n");
    output_c ('\n');
    yyextra->incode = 0;
    yyextra->indent = spacenb (yytext);
    yyextra->type = Octave;
  }
  else
    REJECT;
}

^{WS}*%\}{WS}*$ {
  if (!yyextra->incode) {
    uline();
    // end of Octave documentation comment i.e. "%}"
    output_s ("\n");
    yyextra->ncodes++;
    yyextra->incode = 1;
    yyextra->first = 1;
  }
  else
    REJECT;
}

^#!\/bin\/bash{WS}*$ {
  if (yyextra->line > 2)
    REJECT;
  uline();
}

^:<<'DOC'{WS}*$ {
  if (yyextra->incode) {
    uline();
    // start of Bash documentation comment i.e. ":<<'DOC'"
    if (yyextra->ncodes > 0 && !yyextra->first)
      output_s ("\n~~~\n");
    output_c ('\n');
    yyextra->incode = 0;
    yyextra->indent = spacenb (yytext);
    yyextra->type = Bash;
  }
  else
    REJECT;
}

^DOC{WS}*$ {
  if (!yyextra->incode) {
    uline();
    // end of Bash documentation comment i.e. "DOC"
    yyextra->ncodes++;
    yyextra->incode = 1;
    yyextra->first = 1;
    output_s ("\n");
  }
  else
    REJECT;
}

^{SP}*[\v\n\f] {
  // empty line
  uline();
  if (!yyextra->incode || !yyextra->first)
    output_s (yytext);
}

^{SP}* {
  uline();
  // spaces at the beginning of a line
  if (yyextra->incode) {
    if (yyextra->first && yyextra->type) {
      output_s (yyextra->type);
      output_c ('\n');      
      yyextra->first = 0;
    }
    output_s (yytext);
  }
  else {
    int ns = spacenb (yytext) - yyextra->indent;
    while (ns-- > 0)
      output_c (' ');
  }
}

^{SP}*~~~gnuplot.*$ {
  if (yyextra->incode)
    REJECT;
  uline();
  yyextra->gnuplot = strdup (strstr (yytext, "gnuplot") + 7);
  printf ("<div id=\"plot%d\">\n", yyextra->nplots);
  yyextra->plotype = strdup ("~~~ {.bash}");
  fputs (yyextra->plotype, stdout);
}

^{SP}*~~~pythonplot.*$ {
  if (yyextra->incode)
    REJECT;
  uline();
  yyextra->gnuplot = strdup (strstr (yytext, "pythonplot") + 10);
  printf ("<div id=\"plot%d\">", yyextra->nplots);
  yyextra->plotype = strdup ("~~~ {.python}");
  fputs (yyextra->plotype, stdout);
}

set{SP}+output{SP}*['"][^'"]+['"] |
savefig{SP}*[(]{SP}*['"][^'"]+['"] {
  if (yyextra->gnuplot) {
    uline();
    char * s = strchr (yytext, '\'');
    if (!s)
      s = strchr (yytext, '"');
    s++;
    yyextra->gnuplot_output = strdup (s);
    yyextra->gnuplot_output[strlen(s) - 1] = '\0';
    output_s (yytext);
  }
  else
    REJECT;
}

^{SP}*~~~bib{SP}*$ {
  if (yyextra->inbibtex || yyextra->incode)
    REJECT;
  uline();
  // bibtex file
  char * command = acat ("awk -f $BASILISK/darcsit/hal2bib.awk | "
			 "bibtex2html -a -d -r "
			 "-no-keywords -noabstract -use-keys "
			 "-nodoc -noheader -q | ",
			 "sed -e 's|</table>.*|</table>|' -e '/<\\/table>/q' > ",
			 yyextra->page, ".bib2html",
			 NULL);
  yyextra->inbibtex = popen (command, "w");
  if (yyextra->inbibtex == NULL) {
    perror (command);
    free (command);
    REJECT;
  }
  free (command);
  yyextra->out = yyextra->inbibtex;
}

^{SP}*~~~{SP}*$ {
  uline();
  if (!yyextra->incode && yyextra->gnuplot) {
    if (!yyextra->gnuplot_output) {
      char name[30];
      sprintf (name, "_plot%d.svg", yyextra->nplots);
      yyextra->gnuplot_output = strdup (name);
    }
    output_s ("~~~\n</div>\n");
    output_s ("![");
    char * s = yyextra->gnuplot;
    while (*s != '{' && *s != '\0')
      output_c(*s++);
    printf (" (<a href=\"#\" id=\"buttonplot%d\">script</a>)", yyextra->nplots);
    output_s ("](");
    char timestamp[80];
    sprintf (timestamp, "?%ld", time (NULL));
    char * name = acat (yyextra->basename, "/", yyextra->gnuplot_output, timestamp, NULL);
    output_s (name);
    free (name);
    output_s(")");
    if (*s == '{') {
      while (*s != '}' && *s != '\0')
	output_c(*s++);
      output_c ('}');
    }
    printf ("\n\n<div class=\"plot-script\" id=\"afterplot%d\"></div>", yyextra->nplots);
    free (yyextra->gnuplot);
    yyextra->gnuplot = NULL;
    free (yyextra->gnuplot_output);
    yyextra->gnuplot_output = NULL;
    free (yyextra->plotype);
    yyextra->plotype = NULL;
    yyextra->nplots++;
  }
  else if (yyextra->inbibtex) {
    yyextra->out = stdout;
    pclose (yyextra->inbibtex);
    yyextra->inbibtex = NULL;
    char * file = acat (yyextra->page, ".bib2html", NULL);
    FILE * fp = fopen (file, "r");
    if (!fp)
      perror (file);
    else {
      output_s ("<div class=\"bibtex\">\n");
      int c;
      while ((c = fgetc (fp)) != EOF)
	output_c (c);
      fclose (fp);
      output_s ("</div>\n");
    }
    remove (file);
    free (file);
  }
  else
    output_s (yytext);
}
  
!\[[^\]]*\][(][^)]+\.(png|gif|jpg|mp4|ogv)[)](\([^)]*\))? {
  if (yyextra->incode)
    REJECT;

  uline();
  char * end = strchr (yytext, ']');
  char * name = strchr (end, '(') + 1;

  char * options = strchr (name, ')');
  *options++ = '\0';

  char * link = strdup (name);

  *end = '\0';
  if (!strcmp(link + strlen(link) - 4, ".mp4") ||
      !strcmp(link + strlen(link) - 4, ".ogv")) {
    char * caption = strchr (yytext, '[') + 1,
      * linked = strstr (options, "link");
    if (linked && strchr ("( \t", linked[-1]) && strchr (") \t", linked[4])) {
      // link to video
      for (int i = 0; i < 4; i++)
	linked[i] = ' ';
      output_s ("<div class=\"figure\">");

      output_s ("<a href=\"");
      output_s (link);
      char tstamp[80];
      snprintf (tstamp, 79, "?%ld", time (NULL));
      fputs (tstamp, yyextra->out);
      output_s ("\">");
      
      output_s ("<img ");
      if (options[0] == '(') {
	*strchr (++options, ')') = '\0';
	output_s (options);
	output_c (' ');
      }
      output_s ("src=\"");
      char * snapshot = strdup (link);
      strcpy (snapshot + strlen(snapshot) - 4, ".jpg");
      char * command = malloc (100 + strlen (link) + strlen (snapshot));
      strcpy (command, "ffmpeg -ss 00:00:10 -i '");
      strcat (command, link);
      strcat (command, "' -frames:v 1 -q:v 2 -loglevel quiet -stats -y '");
      strcat (command, snapshot);
      strcat (command, "'");
      int status = system (command); // fixme: warning
      if (status < 0)
	fputs ("", stderr);
      output_s (snapshot);
      output_s ("\">");
      free (snapshot);
      free (command);
      
      output_s ("</a>");
      
      if (caption[0] != '\0') {
	output_s ("<p class=\"caption\">");
	output_s (caption);
	output_s ("</p>");
      }
      output_s ("</div>");
    }
    else { // inline video
      output_s ("<div class=\"figure\"><video ");
      if (options[0] == '(') {
	*strchr (++options, ')') = '\0';
	output_s (options);
	output_c (' ');
      }
      output_s ("controls preload=\"metadata\"><source src=\"");
      output_s (link);
      fprintf (yyextra->out, "?%ld", time (NULL));
      output_s ("\" type = \"video/");
      output_s (!strcmp(link + strlen(link) - 4, ".mp4") ? "mp4" : "ogg");
      output_s ("\"/>Your browser does not support the video tag.</video>");
      if (caption[0] != '\0') {
	output_s ("<p class=\"caption\">");
	output_s (caption);
	output_s ("</p>");
      }
      output_s ("</div>");
    }
  }
  else {
    output_s (yytext);
    output_s ("](");
    output_s (link);
    fprintf (yyextra->out, "?%ld", time (NULL));
    output_c (')');
  }
  free (link);
}

\[[^\[]*\]\({SP}*\) {
  // empty link e.g. [Tutorial]()
  uline();
  char * s = strchr (yytext, ']');
  *s = '\0';
  printf ("[%s](%s)", yytext + 1, yytext + 1);
}

\n {
  uline();
  output_s (yytext);
}

. {
  uline();
  if (yyextra->incode && yyextra->first && yyextra->type) {
    output_s (yyextra->type);
    output_c ('\n');
    yyextra->first = 0;
  }
  output_s (yytext);
}

\"([^\"\\\n]|{ES})*\" {
  /* STRING_LITERAL */
  uline();
  output_s (yytext);
}

%%

static char * append (char * s, char * s1)
{
  if (s) {
    if (strstr (s, s1))
      return s;
    char * n = acat (s, "<br>", s1, NULL);
    free (s);
    return n;
  }
  return strdup (s1);
}

static int scan_errors (FILE * fp, char * file, struct MyScanner * scan,
			int output)
{
  char * header;
  header = acat (file, ".c:", NULL);
  char * line = NULL;
  size_t n = 0, nl = 0;
  while (getline (&line, &n, fp) > 0) {
    if (!strncmp (line, header, strlen (header)) ||
	!strncmp (line + 1, header, strlen (header))) {
      char * lineno = line; while (*lineno != ':') lineno++;
      char * type = ++lineno; while (*type != ':' && *type != '\0') type++;
      if (*type == ':') {
	*type++ = '\0';
	if (*type >= '0' && *type <= '9') {
	  while (*type != ':' && *type != '\0') type++;
	  *type++ = '\0';
	}
	while (strchr(" \t", *type)) type++;
	char * msg = type; while (*msg != ':' && *msg != '\0') msg++;
	if (*msg == ':') {
	  *msg++ = '\0';
	  if (msg[strlen(msg)-1] == '\n')
	    msg[strlen(msg)-1] = '\0';
	  if (!strcmp (type, "error") || !strcmp (type, "fatal error") ||
	      !strcmp(type, "warning")) {
	    int line = atoi (lineno);
	    if (line > 0) {
	      if (line > scan->nerror) {
		scan->error = realloc (scan->error, line*sizeof (Error));
		int i;
		for (i = scan->nerror; i < line; i++)
		  scan->error[i].error = scan->error[i].warning = NULL;
		scan->nerror = line;
	      }
	      if (!strcmp (type, "error") || !strcmp (type, "fatal error"))
		scan->error[line-1].error = 
		  append (scan->error[line-1].error, msg);
	      else
		scan->error[line-1].warning = 
		  append (scan->error[line-1].warning, msg);
	    }
	  }
	  free (line); line = NULL;
	}
	else if (type[0] != '\n')
	  line[strlen(line)] = ':';
	else {
	  free (line); line = NULL;
	}	  
      }
    }
    if (output && line) {
      if (nl < 20) {
	char * s = line;
	while (*s != '\0') {
	  if (*s == '`')
	    *s = '\'';
	  s++;
	}
	if (line[strlen(line)-1] == '\n')
	  line[strlen(line)-1] = '\0';
	output_s1 (scan, line);
	output_s1 (scan, "<br>");
      }
      else if (nl == 20)
	output_s1 (scan, "...<br>");
      nl++;
    }
    free (line); line = NULL;
  }
  free (header);

  return scan->nerror;
}

static int root_is (const char * path, const char * root)
{
  return !strncmp (path, root, strlen(root)) &&
    (path[strlen(root)] == '\0' || path[strlen(root)] == '/');
}

static int tail_is (const char * path, const char * tail)
{
  int plen = strlen (path);
  int tlen = strlen (tail);
  return plen >= tlen && !strcmp (&path[plen - tlen], tail);
}

static void usage (struct MyScanner * scan)
{
  char * s = acat (scan->page, ".itags", NULL);
  FILE * fp = fopen (s, "r");
  if (!fp)
    perror (s);
  else {
    char type[10], title[80], file[80], line[10];
    int nf = 0, ne = 0, nt = 0;
    while (fscanf (fp, "%s %[^\t]\t%s %s", type, title, file, line) == 4) {
      nf++;
      if (root_is (file, "/src/test")) nt++;
      if (root_is (file, "/src/examples")) ne++;
    }
    if (nf > 0) {
      output_s1 (scan, "## Usage\n");
      rewind (fp);
      while (fscanf (fp, "%s %[^\t]\t%s %s", type, title, file, line) == 4)
	if (!root_is (file, "/src/test") &&
	    !root_is (file, "/src/examples")) {
	  output_s1 (scan, "* ["); output_s1 (scan, title); 
	  output_s1 (scan, "](");
	  output_s1 (scan, file); output_s1 (scan, ")\n");
	}
      if (ne > 0) {
	output_s1 (scan, "\n### Examples\n");
	rewind (fp);
	while (fscanf (fp, "%s %[^\t]\t%s %s", type, title, file, line) == 4)
	  if (root_is (file, "/src/examples")) {
	    output_s1 (scan, "* ["); output_s1 (scan, title); 
	    output_s1 (scan, "](");
	    output_s1 (scan, file); output_s1 (scan, ")\n");
	  }
      }
      if (nt > 0) {
	output_s1 (scan, "\n### Tests\n");
	rewind (fp);
	while (fscanf (fp, "%s %[^\t]\t%s %s", type, title, file, line) == 4)
	  if (root_is (file, "/src/test")) {
	    output_s1 (scan, "* ["); output_s1 (scan, title); 
	    output_s1 (scan, "](");
	    output_s1 (scan, file); output_s1 (scan, ")\n");
	  }
      }
    }
    fclose (fp);
  }
  free (s);
}

void literate (FILE * fp, const char * page, int code)
{
  yyscan_t scanner;
  struct MyScanner sdata;
  sdata.in = fp;
  sdata.inbibtex = NULL;

  if (code) {
    if (tail_is (page, ".c") || tail_is (page, ".h"))
      sdata.type = C;
    else if (tail_is (page, ".m"))
      sdata.type = Octave;
    else if (tail_is (page, ".py"))
      sdata.type = Python;
    else
      sdata.type = Bash;
    sdata.first = 1;
    sdata.ncodes = 1;
  }
  else {
    sdata.type = NULL;
    sdata.first = 0;
    sdata.ncodes = 0;
  }
  sdata.incode = code;

  sdata.plotype = NULL;
  sdata.indent = 0;
  sdata.out = stdout;
  sdata.i = 0;
  sdata.error = NULL;
  sdata.nerror = sdata.nplots = 0;
  sdata.line = 1;
  sdata.page = strdup (page);
  sdata.basename = strdup (page);  
  char * s = &(sdata.basename[strlen(sdata.basename) - 1]);
  while (*s != '.' && s != sdata.basename) s--;
  if (*s == '.')
    *s = '\0';
  sdata.gnuplot = sdata.gnuplot_output = NULL;

  if (tail_is (page, ".c")) {
    char * page1 = strdup (page);
    char * file = &page1[strlen(page1) - 2];
    *file-- = '\0';
    while (file != page1 && *file != '/')
      file--;
    if (file != page1)
      *file++ = '\0';

    // scan for errors
    {
      char * fail = acat (file, "/fail", NULL);
      FILE * fp = fopen (fail, "r");
      if (fp) {
	scan_errors (fp, file, &sdata, 0);
	fclose (fp);
      }
      free (fail);
    }

    // scan for warnings
    {
      char * warn = acat (file, "/warn", NULL);
      FILE * fp = fopen (warn, "r");
      if (fp) {
	scan_errors (fp, file, &sdata, 0);
	fclose (fp);
      }
      free (warn);
    }
    
    free (page1);
  }

  yylex_init_extra (&sdata, &scanner);
  yyset_out (stdout, scanner);
  yyset_debug (2, scanner);
  yylex (scanner);
  yylex_destroy (scanner);

  if (sdata.ncodes > 0 && sdata.incode && !sdata.first)
    output_s1 ((&sdata), "\n~~~\n");

  if (sdata.ncodes > 0 && tail_is (page, ".h"))
    usage (&sdata);

#if 0
  fputc ('"', stderr);
  fputs (sdata.output, stderr);
  fputs ("\"\n", stderr);
#endif

  int i;
  for (i = 0; i < sdata.nerror; i++) {
    free (sdata.error[i].error);
    free (sdata.error[i].warning);
  }
  free (sdata.error);
  free (sdata.page);
  free (sdata.basename);
}

int main (int argc, char * argv[])
{
  if (argc != 3) {
    fprintf (stderr, "usage: ./literate FILE CODE\n");
    return 1;
  }
  char * name = acat (argv[1], ".page", NULL);
  FILE * f = fopen (name, "r");
  if (!f)
    f = fopen (argv[1], "r");
  if (!f) {
    perror (name);
    return 1;
  }
  free (name);

  literate (f, argv[1], atoi(argv[2]));
  return 0;
}

