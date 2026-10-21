/**
# An embedded gview

This is a version of [gview](gview/display.c) which can be embedded
directly in a simulation. It replaces the client/server model with a
simple set of two queues within the same process: this transfer
mechanism is implemented in [gview/transfer_local.c]().

For the moment this has little advantage over using the standard
client/server system (using WebSockets), but this will be useful when
we implement direct GPU rendering, which require the client and server
to run on the same process. */

#include "gview/display.h"
#include <pthread.h>

static pthread_t gview_thread;

static void * gview_main (void *arg)
{
  // display_set_debug (true);
  if (display_main ("local", 0, NULL))
    exit (1);
  return NULL;
}

static void gview_terminate()
{
  if (pid() == 0) {
    display_quit();
    pthread_join (gview_thread, NULL);
  }
}

@init_solver
void gview_init()
{
  if (pid() == 0 &&
      pthread_create (&gview_thread, NULL, gview_main, NULL) != 0) {
    perror ("pthread_create");
      exit (1);
  }
  free_solver_func_add (gview_terminate);
}

#define DISPLAY_JS ""
#pragma autolink $BASILISK/gview/transfer_local.o -L$BASILISK/gview -ldisplay \
  -Lgl -lglutils -lglfw -limgui -lstb -lstdc++ -lm
#define DISPLAY -1
#include "display.h"
