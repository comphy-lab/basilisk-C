%option noyywrap yylineno
%{

  static int is_allowed (const char * key) {
    // White list of HTML tags
    static char * allowed[20] = {
      "a", "br", "caption", "center", "div", "hr", "img", "p", "pre",
      "source", "span", "table", "td", "tr", "video",
      NULL
    };
    char ** i = allowed;
    while (*i) {
      if (!strcmp (key, *i))
	return 1;
      i++;
    }
    return 0;
  }
  
  int incode = 0, inmath = 0, inscript = 0;

  char * filename = "";
  
  static void echo() {
    fputs (yytext, yyout);
  }

%}

ID     [a-zA-Z0-9_]
SP     [ \t]
WS     [ \t\v\n\f]
ES     (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))

%%

\"([^\"\\\n]|{ES})*\"  echo(); /* STRING_LITERAL */

{WS}*\/[*][*]{SP}* {
  // start of C documentation comment i.e. "/**"
  incode = 0;
  echo();
}

{SP}*[*]\/{SP}* {
  // end of any C comment block i.e. "*/"
  if (!incode)
    incode = 1;
  echo();
}

{WS}*\"\"\"{SP}* {
  if (incode)
    // start of Python documentation comment i.e. """
    incode = 0;
  else
    // end of Python documentation comment
    incode = 1;
  echo();
}

^{WS}*%\{{WS}*$ {
  if (incode) {
    // start of Octave documentation comment i.e. "%{"
    incode = 0;
    echo();
  }
  else
    REJECT;
}

^{WS}*%\}{WS}*$ {
  if (!incode) {
    // end of Octave documentation comment i.e. "%}"
    incode = 1;
    echo();
  }
  else
    REJECT;
}

^#!\/bin\/bash{WS}*$ {
  if (yylineno > 2)
    REJECT;
}

^:<<'DOC'{WS}*$ {
  if (incode) {
    // start of Bash documentation comment i.e. ":<<'DOC'"
    incode = 0;
    echo();
  }
  else
    REJECT;
}

^DOC{WS}*$ {
  if (!incode) {
    // end of Bash documentation comment i.e. "DOC"
    incode = 1;
    echo();
  }
  else
    REJECT;
}

[<][/]{0,1}[a-zA-Z]+({SP}[^<>$]+){0,1}[>] {
  if (!incode && !inmath && !inscript) {
    char * key = strdup (yytext + 1), * s = key;
    while (!strchr(" \t>", *s)) s++;
    *s = '\0';
    // fprintf (stderr, "%s:%d: %s\n", filename, yylineno, yytext);
    if (!strcmp(key, "pre"))
      inscript = 1;
    char * skey = key[0] == '/' ? key + 1 : key;
    if (is_allowed (skey))
      echo();
    else
      fprintf (stderr, "%s:%d: HTML tag '%s' is not allowed\n",
	       filename, yylineno, skey);
    free (key);
  }
  else
    REJECT;
}

[<]/pre[>] {
  if (inscript)
    inscript = 0;
  else
    REJECT;
}

[$] {
  if (incode)
    REJECT;
  inmath = !inmath;
  echo();
}

^[~]{3} {
  if (incode)
    REJECT;
  inscript = !inscript;
  echo();
}

\'.\' {
  echo(); // quoted character
}

.                   echo();
[\n]                echo();

%%

int main (int argc, char * argv[])
{
  if (argc > 1)
    filename = argv[1];
  return yylex();
}
