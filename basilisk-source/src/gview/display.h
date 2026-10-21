#ifndef DISPLAY_H
#define DISPLAY_H

#include <stdbool.h>

void display_add_item (const char * command);
void display_remove_item (unsigned int id);
void display_update_item (unsigned int id, const char * command);
void display_item_visible (unsigned int id, bool visible);
void display_connect (const char * host, int port);
void display_disconnect();
void display_quit();
bool display_load_data (const char *filename, const char *data);
bool display_load (const char * filename);

void display_reset_view ();
void display_front_view ();
void display_back_view ();
void display_left_view ();
void display_right_view ();
void display_top_view ();
void display_bottom_view ();
void display_iso_view ();

void display_set_running (int play);
void display_set_connected (int connected);
void display_set_debug (bool debug);
int display_main (const char * host, int port, const char * gv);

typedef struct {
  const char * name;
  char ** keys;
  char ** values;
  int n;
} DisplayControlStrings;

void display_control_command (const char * command);

typedef struct {
  float tx, ty, tz;
  float sx, sy, sz;
  float rotation[4];       /* x, y, z, w, as in view.h */
  float fov;
  float near, far;
  float bg[4];
  
  double mouse_x, mouse_y;
  bool rotating;
  bool panning;
  bool zooming;
} DisplayCamera;

extern DisplayCamera Camera;

#endif
