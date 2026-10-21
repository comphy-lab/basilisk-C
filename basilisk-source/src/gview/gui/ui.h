#include <GLFW/glfw3.h>
#include "../display.h"

void ui_init (GLFWwindow *window);
bool ui_add_item (const char * command);
bool ui_display (GLFWwindow *window);
void ui_viewport (GLFWwindow *window, int * xpos, int * ypos, int * width, int * height);
void ui_set_window_size (GLFWwindow *window, int width, int height);

bool ui_want_capture_mouse();
bool ui_want_capture_keyboard();

void ui_set_running (int running);
void ui_set_connected (int connected);
void ui_set_dimension (int dim);
void ui_set_host (const char * host, int port);
void ui_set_error (int id, const char * error);
void ui_copy_error (int id, int src);
void ui_set_status (const char * status);
  
void ui_add_control (DisplayControlStrings control);
void ui_destroy();
