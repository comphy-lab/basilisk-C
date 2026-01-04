// Returns successfully if a file starts with /** or %{ or """ or :<<'DOC'

#include <stdio.h>
#include <string.h>

static int pattern (FILE * fp, char * pat)
{
  int c = fgetc (fp);
  do {
    while (c != EOF && strchr(" \t", c))
      c = fgetc (fp);
    char * s = pat;
    while (*s != '\0' && c == *s++)
      c = fgetc (fp);
    if (*s == '\0') {
      while (c != EOF && strchr(" \t", c))
	c = fgetc (fp);
      if (strchr("\v\n\f\r", c))
	return 0;      
    }
    while (c != EOF && !strchr("\v\n\f\r", c))
      c = fgetc (fp);
    c = fgetc (fp);
  } while (c != EOF);
  return 1;
}
    
int main (int argc, char * argv[])
{
  if (argc != 2) {
    fprintf (stderr, "usage: pagemagic FILE\n");
    return 1;
  }
  FILE * fp = fopen (argv[1], "r");
  if (fp == NULL) {
    perror (argv[1]);
    return 1;
  }
  int len = strlen(argv[1]);
  if (len > 2 && (!strcmp (argv[1] + len - 2, ".c") ||
		  !strcmp (argv[1] + len - 2, ".h")))
    return pattern (fp, "/**");
  else if (len > 2 && !strcmp (argv[1] + len - 2, ".m"))
    return pattern (fp, "%{");
  else if (len > 3 && !strcmp (argv[1] + len - 3, ".py"))
    return pattern (fp, "\"\"\"");
  else
    return pattern (fp, ":<<'DOC'");
  return 1;
}
