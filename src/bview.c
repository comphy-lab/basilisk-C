/**
# Basilisk View "dump files" server

This can be used to interactively display [dump files](output.h#dump)
using the [javascript interface](jview/README). 

Typical usage is:

~~~bash
bview3D /path/to/file/dump
~~~

(note that the appropriate 2D/3D version must be chosen). The
connection URL for the client will then be displayed. */

#define DISPLAY 0
#define DISPLAY_NO_CONTROLS
#include "display.h"

int main (int argc, char * argv[])
{
  char * file = argc > 1 ? argv[1] : "dump";

  if (!restore (file = file, list = all)) {
    fprintf (stderr, "%s: could not restore from '%s'\n", argv[0], file);
    exit (1);
  }
  restriction (all);
  fields_stats();

  fputc ('\n', stderr);
  display_url (stderr);
  fputc ('\n', stderr);

  display ("box();");
  while (1) {
    if (display_poll (-1))
      display_update (INT_MAX);
  }
  display_destroy();
}
