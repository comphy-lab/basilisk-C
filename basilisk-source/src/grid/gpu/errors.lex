%option noyywrap
%option yylineno
%option noinput
%{

  /**
  # Error message parsing for compilation of GLSL code on GPUs
  */
  static char * sout;
  static char * fname;
  static int line;
  
  static char * str_append (char * dst, const char * src) {
    if (dst == NULL) {
      dst = malloc (strlen(src) + 1);
      dst[0] = '\0';
    }
    else
      dst = realloc (dst, strlen(dst) + strlen(src) + 1);
    strcat (dst, src);
    return dst;
  }

  static int atoi_check (char * start, char * end)
  {
    char * s = start;
    while (s != end) {
      if (*s < '0' || *s > '9')
	return -1;
      s++;
    }
    *end = '\0';
    return atoi (start);
  }

  /* parse Intel GPU error messages of the form:
     0:31(15):... */
  static int parse_intel (const char * errors, const char ** msg)
  {
    char * end = strchr (errors, '\n');
    char * sline = strchr (errors, ':');
    if (!sline) return -1; sline++;
    if (end && end < sline) return -1;
    char * s1 = strchr (sline, '(');
    if (!s1) return -1;
    *msg = strchr (s1, ':');
    if (!*msg) return -1; (*msg)++;
    return atoi_check (sline, s1);
  }

  /* parse Nvidia GPU error messages of the form:
     0(31) :... */
  static int parse_nvidia (const char * errors, const char ** msg)
  {
    char * end = strchr (errors, '\n');
    char * sline = strchr (errors, '(');  
    if (!sline) return -1; sline++;
    if (end && end < sline) return -1;
    char * s1 = strchr (sline, ')');
    if (!s1) return -1;
    *msg = strchr (s1, ':');
    if (!*msg) return -1; (*msg)++;
    return atoi_check (sline, s1);
  }

  static char * next_error (const char * errors, int * line, const char ** msg)
  {
    if (!errors || *errors == '\0')
      return NULL;
    
    *line = parse_intel (errors, msg);
    if (*line < 0)
      *line = parse_nvidia (errors, msg);
    if (*line < 0) {
      // unknown error format: assumes first line
      *line = 1;
      *msg = errors;
    }

    while (strchr (" \t", **msg)) (*msg)++;
    
    char * s3 = strchr (*msg, '\n');
    if (s3) {
      *s3 = '\0';
      s3++;
    }

    return s3;
  }
  
  static char * merge (const char * source, const char * errors)
  {
    char * merged = NULL;
    const char * msg;
    int eline, line = 1;
    char * error = next_error (errors, &eline, &msg);
    char * src = strdup (source), * s = src;
    while (*s != '\0') {
      if (*s == '\n') {
	while (error && eline <= line) {
	  merged = str_append (merged, "@");
	  if (eline != line)
	    merged = str_append (merged, "!");
	  merged = str_append (merged, msg);
	  merged = str_append (merged, "@");
	  error = next_error (error, &eline, &msg);
	}
	line++;
      }
      char c[] = " ";
      c[0] = *s;
      merged = str_append (merged, c);
      s++;
    }
    free (src);
    return merged;
  }
%}

ID     [a-zA-Z0-9_]
SP     [ \t]
WS     [ \t\v\n\f]
ES     (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))

%%

^"//"{SP}*#{SP}*(line)?{SP}+[0-9]+({SP}+\".*\")? {
  char * sline = yytext;
  while (!strchr ("0123456789", *sline)) sline++;
  line = atoi (sline) - 1;

  if (strchr (yytext, '"')) {
    fname = sline; while (*fname >= '0' && *fname <= '9') fname++;
    *fname = '\0'; fname++;
    fname = strchr (fname, '"'); fname++;
    *strchr (fname, '"') = '\0';
  }
}

@[^@]*@ {
  if (0) // !strncmp (yytext, "@error ", 7))
    yytext += 6;
  else {
    if (fname)
      sout = str_append (sout, fname);
    else
      sout = str_append (sout, "(fragment shader)");
    sout = str_append (sout, ":");
    char s[81];
    snprintf (s, 80, "%d", line);
    sout = str_append (sout, s);
    sout = str_append (sout, ": GLSL: ");
  }
  *strchr (yytext + 1, '@') = '\0';
  sout = str_append (sout, yytext + 1);
  sout = str_append (sout, "\n");
}

\n { line++; }

. ; // very important!!

%%

char * gpu_errors (const char * errors, const char * source, char * fout)
{
  if (0) yyunput (0, NULL); // just prevents 'yyunput unused' compiler warning

  char * merged = merge (source, errors);
  sout = fout;
  sout = str_append (NULL, "");
  YY_BUFFER_STATE buf = yy_scan_string (merged);
  fname = NULL;
  line = 1;
  yylex();
  yy_delete_buffer (buf);
  free (merged);

  return sout;
}
