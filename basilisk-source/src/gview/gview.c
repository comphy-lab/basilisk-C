#include "display.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void usage (const char * command)
{
  fprintf (stderr,
           "usage: %s [OPTIONS] [HOST] [PORT]\n"
           "   The default HOST is 'localhost' and default PORT is 7100.\n"
           "   The OPTIONS are:\n"
           "   -g FILE     load .gv parameters\n"
           "   -d, --debug debugging output\n"
           "   -h, --help  this message\n",
           command);
}

int main (int argc, char * argv[])
{
  int port = 7100, hosti = -1, porti = -1;
  const char * host = "localhost", * gv = NULL;
  for (int i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "--debug") || !strcmp (argv[i], "-d"))
      display_set_debug (true);
    else if (!strcmp (argv[i], "--help") || !strcmp (argv[i], "-h")) {
      usage (argv[0]);
      return 0;
    }
    else if (!strcmp (argv[i], "-g")) {
      if (i == argc - 1) {
        fprintf (stderr, "%s: error: missing FILE\n", argv[0]);
        usage (argv[0]);
        return 1;
      }
      gv = argv[i+1];
      i++;
    }
    else if (argv[i][0] == '-') {
      fprintf (stderr,
               "%s: unknown option '%s'\n"
               "  try %s -h\n",
               argv[0], argv[i], argv[0]);
      return 1;
    }
    else if (hosti < 0)
      hosti = i;
    else if (porti < 0)
      porti = i;
    else {
      usage (argv[0]);
      return 1;
    }
  }
  if (hosti > 0)
    host = argv[hosti];
  if (porti > 0)
    port = atoi (argv[porti]);

  return display_main (host, port, gv);
}
