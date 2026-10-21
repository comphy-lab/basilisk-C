/**
# A GUI for gview

## To Do

* Disconnect button doesn't do anything yet
* Error messages are not displayed
*/

#include "bindings/dcimgui.h"
#include "bindings/dcimgui_impl_glfw.h"
#include "bindings/dcimgui_impl_opengl3.h"
#include "ui.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef __EMSCRIPTEN__
# include <emscripten.h>
#endif

/**
# Configuration */

static const float font_size = 18.f;
static const float ui_top = 21.f;

/**
# GUI state */

typedef struct {
  char * name, * tooltip, * type, * value, * min, * max;
  union {
    int i;
    float f;
  } val;
} UiControl;

static struct {
  GLFWwindow *window;
  UiControl * controls;
  float draw_window_width;
  int running;   // -1: paused   0: running       1: step       -2: no controls
  int connected; // -1: waiting  0: disconnected  1: connected
  int dimension;
  char host[80], * status;
  int port;
} UIState = {
  .draw_window_width = 440,
  .host = "localhost", .port = 7100,
  .dimension = 2
};

void ui_set_host (const char * host, int port)
{
  strncpy (UIState.host, host, 79);
  UIState.port = port;
}

void ui_set_status (const char * status)
{
  free (UIState.status);
  UIState.status = status ? strdup (status) : NULL;
}

/**
# Draw item */

typedef struct UIDrawItem UIDrawItem;

struct UIDrawItem {
  unsigned int id;
  const char *name;
  char * error;
  bool visible;

  bool (*properties) (UIDrawItem *);
  const char *(*summary) (UIDrawItem *);
  char *(*command) (UIDrawItem *);
};

/**
# Hyperlinks */

static bool imgui_link (const char *label, const char *url)
{
  static const ImVec4 normal = {0.2f, 0.5f, 1.f, 1.f};
  static const ImVec4 hover = {0.4f, 0.7f, 1.f, 1.f};

  ImVec2 pos = IgGetCursorScreenPos ();
  ImVec2 size = IgCalcTextSize (label);

  ImVec2 min = pos;
  ImVec2 max = {pos.x + size.x, pos.y + size.y};

  bool hovered = IgIsMouseHoveringRect (min, max);

  ImU32 color = IgColorConvertFloat4ToU32 (hovered ? hover : normal);

  IgPushStyleColor (ImGuiCol_Text, color);
  IgTextUnformatted (label);
  IgPopStyleColor ();

  if (hovered) {
    IgSetMouseCursor (ImGuiMouseCursor_Hand);
  }

  ImDrawList *draw_list = IgGetWindowDrawList ();

  ImVec2 p1 = {min.x, max.y};
  ImVec2 p2 = {max.x, max.y};

  ImDrawList_AddLine (draw_list, p1, p2, color);

  if (hovered && IgIsMouseClicked (ImGuiMouseButton_Left)) {
#ifdef __EMSCRIPTEN__
    EM_ASM ({ window.open (UTF8ToString ($0), "_blank"); }, url);
#else
    char command[1024];

    snprintf (command, sizeof command, "xdg-open '%s' &", url);
    if (system (command))
      fprintf (stderr, "%s:%d: %s failed\n", __FILE__, __LINE__, command);
#endif
    return true;
  }

  return false;
}

/**
## Automated interface definition

We first add a few definitions which need to match those in
[/src/draw.h](). */

typedef int Colormap;

static const Colormap
  jet = 0,
  cool_warm = 1,
  gray = 2,
  randomap = 3,
  blue_white_red = 4;

static const int ncolormaps = 5;
static const char * colormaps[5] = {
  "jet", "cool_warm", "gray", "randomap", "blue_white_red"
};

typedef struct { double x, y, z; } coord;

typedef void * MapFunc, * bview;

#define HUGE 1e30

/**
These files are generated automatically by [/src/Makefile](). */

#include "../../draw_ui.h"
#include "../../parse.h"

/**
# Controls */

static
void new_control (DisplayControlStrings * control, UiControl * c)
{
  memset (c, 0, sizeof (UiControl));
  c->name = strdup (control->name);
  for (int i = 0; i < control->n; i++)
    if (!strcmp (control->keys[i], "tooltip")) {
      int np = 0, len = 0;
      for (char * s = control->values[i]; *s != '\0'; s++, len++)
        if (*s == '%')
          np++;
      if (!np)
        c->tooltip = strdup (control->values[i]);
      else {
        char * s1 = malloc (len + np + 1);
        c->tooltip = s1;
        for (char * s = control->values[i]; *s != '\0'; s++, s1++)
          if (*s == '%') *s1++ = '%', *s1 = '%';
          else *s1 = *s;
        *s1 = '\0';
      }
    }
    else if (!strcmp (control->keys[i], "type"))
      c->type = strdup (control->values[i]);
    else if (!strcmp (control->keys[i], "value"))
      c->value = strdup (control->values[i]);
    else if (!strcmp (control->keys[i], "min"))
      c->min = strdup (control->values[i]);
    else if (!strcmp (control->keys[i], "max"))
      c->max = strdup (control->values[i]);
  if (c->type && c->value) {
    if (!strcmp (c->type, "int"))
      c->val.i = atoi(c->value);
    else if (!strcmp (c->type, "double"))
      c->val.f = atof(c->value);
  }
}

static void update_control (const UiControl * c)
{
  assert (c->type);
  char * s = NULL;
  if (!strcmp (c->type, "int"))
    assert (asprintf (&s, "#%s:%d", c->name, c->val.i) > 0);
  else if (!strcmp (c->type, "double"))
    assert (asprintf (&s, "#%s:%g", c->name, c->val.f) > 0);
  else
    assert (false);
  display_control_command (s);
  free (s);
}

static void destroy_control (UiControl * c)
{
  free (c->name);
  free (c->tooltip);
  free (c->type);
  free (c->value);
  free (c->min);
  free (c->max);
}

static void destroy_controls()
{
  for (UiControl * c = UIState.controls; c && c->name; c++)
    destroy_control (c);
  free (UIState.controls);
  UIState.controls = NULL;
}

void ui_add_control (DisplayControlStrings control)
{
  int n = 0;
  for (UiControl * c = UIState.controls; c && c->name; c++, n++)
    if (!strcmp (c->name, control.name)) {
      const char * type = NULL, * value = NULL;
      for (int i = 0; i < control.n; i++)
        if (!strcmp (control.keys[i], "type"))
          type = control.values[i];
        else if (!strcmp (control.keys[i], "value"))
          value = control.values[i];
      if (type && value && !strcmp (c->type, type)) {
        if (!strcmp (c->type, "int"))
          c->val.i = atoi(value);
        else if (!strcmp (c->type, "double"))
          c->val.f = atof(value);
        return;
      }
    }
  UIState.controls = realloc (UIState.controls, (n + 2)*sizeof (UiControl));
  UIState.controls[n + 1].name = NULL;
  new_control (&control, &UIState.controls[n]);
}

/**
# Draw list items */

typedef struct {
  UIDrawItem **items;
  int nitems;
} UIDrawList;

static UIDrawList draw_list = {0};
static UIDrawItem * selected_draw_item = NULL, * view_item = NULL;
static bool properties_open;

static void ui_add_draw_item (UIDrawItem *item)
{
  draw_list.items = realloc (draw_list.items,
                             (draw_list.nitems + 1)*
                             sizeof(UIDrawItem *));
  item->id = draw_list.nitems - 1;
  draw_list.items[draw_list.nitems++] = item;
  if (item != view_item) {
    char * command = item->command (item);
    display_add_item (command);
    free (command);
  }
}

void ui_set_error (int id, const char * error)
{
  assert (id <= draw_list.nitems);
  free (draw_list.items[id + 1]->error);
  if (error)
    draw_list.items[id + 1]->error = strdup (error);
  else
    draw_list.items[id + 1]->error = NULL;
}

void ui_copy_error (int id, int src)
{
  assert (src <= draw_list.nitems);
  ui_set_error (id, draw_list.items[src + 1]->error);
}

static void free_item (UIDrawItem *item)
{
  display_remove_item (item->id);
  free (item->error);
  free (item);
}

static void ui_delete_draw_item (UIDrawItem * item)
{
  if (selected_draw_item == item) {
    selected_draw_item = NULL;
    properties_open = false;
  }

  unsigned int index;
  for (index = 0; index < draw_list.nitems; index++)
    if (draw_list.items[index] == item)
      break;
  assert (index < draw_list.nitems);
  
  free_item (draw_list.items[index]);

  for (unsigned int i = index; i < draw_list.nitems - 1; i++) {
    draw_list.items[i] = draw_list.items[i + 1];
    draw_list.items[i]->id = i - 1;
  }
  draw_list.nitems--;

  if (draw_list.nitems == 0) {
    free (draw_list.items);
    draw_list.items = NULL;
  }
}

/**
# Save commands */

#ifdef __EMSCRIPTEN__
static void browser_download (const char *filename, const char *data)
{
  EM_ASM ({
    const filename = UTF8ToString ($0);
    const data = UTF8ToString ($1);

    const blob = new Blob ([data], { type: "text/plain" });
    const url = URL.createObjectURL (blob);

    const a = document.createElement ("a");
    a.href = url;
    a.download = filename;
    a.click ();

    URL.revokeObjectURL (url);
  }, filename, data);
}

EMSCRIPTEN_KEEPALIVE
void gview_file_loaded (const char *filename, const char *data)
{
  display_load_data (filename, data);
}

#else
# include "tinyfiledialogs.h"
#endif

static
const char *ui_commands (void)
{
  static char *commands;
  static size_t size;

  size_t length = 0;

  if (!commands) {
    size = 1024;
    commands = malloc (size);
  }

  commands[0] = '\0';
  for (int i = 0; i < draw_list.nitems; i++) {
    UIDrawItem *item = draw_list.items[i];
    if (!item->visible)
      continue;
    char *command = item->command (item);
    if (!command)
      continue;
    size_t n = strlen (command);
    if (length + n + 2 > size) {
      while (length + n + 2 > size)
        size *= 2;
      commands = realloc (commands, size);
    }
    memcpy (&commands[length], command, n);
    length += n;
    commands[length++] = '\n';
    commands[length] = '\0';
    free (command);
  }

  return commands;
}

#ifdef __EMSCRIPTEN__
static void gview_save (const char *filename)
{
  const char *commands = ui_commands ();
  if (commands)
    browser_download (filename, commands);
}
#else
static void gview_save (const char *filename)
{
  const char *commands = ui_commands();
  if (!commands)
    return;
  FILE *fp = fopen (filename, "w");
  if (!fp) {
    perror (filename);
    return;
  }
  fputs (commands, fp);
  fclose (fp);
}
#endif

/**
# Draw list */

typedef struct {
  const char * name;
  UIDrawItem * (* new)();
} DrawItemMenu;

void ui_set_dimension (int dim)
{
  UIState.dimension = dim;
}

static bool ui_draw_menu (void)
{
  bool changed = false;
  if (IgBeginMenu ("Draw")) {
    static const DrawItemMenu items[] = {
      { "squares",  ui_squares_new },
      { "draw_vof", ui_draw_vof_new },
      { "cells",    ui_cells_new },
      { "box",      ui_box_new },
      { NULL }
    };
    /* 2D only */
    static const DrawItemMenu items2D[] = {
      { "isoline",  ui_isoline_new },
      { "labels",   ui_labels_new },
      { "vectors",  ui_vectors_new },
      { "lines",    ui_lines_new },      
      { NULL }
    };
    /* 3D only */
    static const DrawItemMenu items3D[] = {
      { "isosurface", ui_isosurface_new },
      { NULL }
    };
    for (const DrawItemMenu * i = items; i->name; i++)
      if (IgMenuItem (i->name)) {
        UIDrawItem * item = i->new();
        ui_add_draw_item (item);
        selected_draw_item = item;
        properties_open = true;
        changed = true;
      }
    if (UIState.dimension < 3) {
      for (const DrawItemMenu * i = items2D; i->name; i++)
        if (IgMenuItem (i->name)) {
          UIDrawItem * item = i->new();
          ui_add_draw_item (item);
          selected_draw_item = item;
          properties_open = true;
          changed = true;
        }
    }
    else {
      for (const DrawItemMenu * i = items3D; i->name; i++)
        if (IgMenuItem (i->name)) {
          UIDrawItem * item = i->new();
          ui_add_draw_item (item);
          selected_draw_item = item;
          properties_open = true;
          changed = true;
        }
    }
    IgEndMenu ();
  }
  return changed;
}

static void ui_connect()
{
  assert (!UIState.connected);
  display_connect (UIState.host, UIState.port);
}

/**
# Properties */

static bool ui_draw_properties (void)
{
  if (!properties_open || !selected_draw_item)
    return false;

  IgBeginDisabled (UIState.connected != 1);
  
  bool changed = false;

  ImVec2 avail = IgGetContentRegionAvail();
  if (UIState.status)
    avail.y -= IgGetTextLineHeight() + 4;
  
  if (IgBeginChild (selected_draw_item->name, avail, true, 0)) {
    if (selected_draw_item->error) {
      IgPushStyleColor (ImGuiCol_Text, IM_COL32(255,0,0,255));
      IgTextUnformatted (selected_draw_item->error);
      IgPopStyleColor();
    }
    changed = selected_draw_item->properties (selected_draw_item);
    if (changed && selected_draw_item != view_item) {
      char * command = selected_draw_item->command (selected_draw_item);
      display_update_item (selected_draw_item->id, command);
      free (command);
    }
  }

  IgEndChild ();

  IgEndDisabled();

  return changed;
}

static bool ui_draw_list (int width, int height)
{
  bool changed = false;

  static float draw_window_width_start;
  static bool draw_window_width_initialized = false;
  int lwidth = UIState.draw_window_width;

  /**
  # Draw window */
  
  if (!draw_window_width_initialized &&
      access (IgGetIO ()->IniFilename, R_OK) == 0) {
    /**
    Don't force a size here. Let ImGui restore it from .gview.ini. */
  }
  else {
    IgSetNextWindowPos ((ImVec2){0, ui_top}, ImGuiCond_Always);
    IgSetNextWindowSize ((ImVec2){lwidth, height - ui_top}, ImGuiCond_Always);
  }

  if (IgBegin ("Draw", NULL,
               ImGuiWindowFlags_NoMove |
               ImGuiWindowFlags_NoTitleBar |
               ImGuiWindowFlags_NoResize |
               ImGuiWindowFlags_NoCollapse |
               ImGuiWindowFlags_NoBringToFrontOnFocus)) {

    if (!draw_window_width_initialized) {
      UIState.draw_window_width = IgGetWindowWidth();
      draw_window_width_initialized = true;
    }
    
    /**
    Connection */

    if (strcmp (UIState.host, "local") || UIState.port != 0) {
      if (UIState.connected) {
        if (IgButton ("Disconnect"))
          display_disconnect();
        IgPushStyleColor (ImGuiCol_Text,
                          UIState.connected > 0 ?
                          IM_COL32(0,255,0,255) : // Connected: Green
                          IM_COL32(230,104,30,255)); // Waiting: Orange
      }
      else {
        if (IgButton ("Connect"))
          ui_connect();
        IgPushStyleColor (ImGuiCol_Text, IM_COL32(255,0,0,255)); // Red
      }
      IgSameLine();

#ifdef __EMSCRIPTEN__
      int host_flags = ImGuiInputTextFlags_ReadOnly; // hostname is not editable in browser mode
#else
      int host_flags = UIState.connected ? ImGuiInputTextFlags_ReadOnly : 0;
#endif
      int port_flags = UIState.connected ? ImGuiInputTextFlags_ReadOnly : 0;

      float host_width = IgCalcTextSizeEx (UIState.host, NULL, false, -1.0f).x +
        2.f*IgGetStyle()->FramePadding.x + 4.f;

      if (host_width < 80.f)
        host_width = 80.f;
      
      IgSetNextItemWidth (host_width);
      IgInputText ("##host", UIState.host, sizeof(UIState.host), host_flags);
      
      IgSameLine();
      IgSetNextItemWidth (50.f);
      IgInputScalarEx ("##port", ImGuiDataType_S32, &UIState.port,
                       NULL, NULL, "%d", port_flags);
      
      IgPopStyleColor();
    }

    /**
    ## Display controls */

    if ((UIState.controls ||
         (UIState.connected > 0 && UIState.running > -2))) {
      IgSeparator();
      if (IgBeginTable ("controls", 2,
                        ImGuiTableFlags_BordersInnerV |
                        ImGuiTableFlags_Resizable |
                        ImGuiTableFlags_SizingFixedFit)) {
          
        IgTableSetupColumn ("Name", 0);
        IgTableSetupColumnEx ("Value", ImGuiTableColumnFlags_WidthStretch, 0.f, 0);
      
        if (UIState.connected > 0 && UIState.running > -2) {
          IgTableNextRow();
          IgTableNextColumn ();
          IgTextUnformatted ("Run/Pause");
          IgTableNextColumn ();        
          if (UIState.running != 0) {
            if (IgButton ("RUN"))
              display_set_running (0);
            IgSameLine();
            if (IgButton ("STEP"))
              display_set_running (1);
          }
          else if (IgButton ("PAUSE"))
            display_set_running (-1);
        }

        int i = 1;
        for (UiControl * c = UIState.controls; c && c->name; c++) {
          IgPushIDInt (i++);
          IgTableNextRow();
          IgTableNextColumn();
          IgTextUnformatted (c->name);
          if (c->tooltip && IgIsItemHovered (0))
            IgSetTooltip (c->tooltip);
          if (c->type) {
            bool changed = false;
            if (!strcmp (c->type, "int")) {
              IgTableNextColumn();
              if (IgInputInt ("##value", &c->val.i)) {
                if (c->max && c->val.i > atoi(c->max)) c->val.i = atoi(c->max);
                if (c->min && c->val.i < atoi(c->min)) c->val.i = atoi(c->min);
                changed = true;
              }
            }
            else if (!strcmp (c->type, "double")) {
              IgTableNextColumn();
              if (c->min && c->max)
                IgDragFloatEx ("##value", &c->val.f, 1., atof (c->min), atof(c->max), "%.3f", 0);
              else
                IgInputFloat ("##value", &c->val.f);
              if (IgIsItemDeactivatedAfterEdit())
                changed = true;
            }
            if (changed)
              update_control (c);
          }
          IgPopID ();
        }

        IgEndTable();
      } // Controls table
    }
      
    /**
    ## Draw items */

    IgSeparator();
    
    if (selected_draw_item && selected_draw_item != view_item &&
        IgIsKeyPressed (ImGuiKey_Delete)) {
      ui_delete_draw_item (selected_draw_item);
      changed = true;
    }

    if (IgBeginTable ("draw_items", 3,
                      ImGuiTableFlags_BordersInnerV |
                      ImGuiTableFlags_Resizable |
                      ImGuiTableFlags_SizingFixedFit)) {

      IgTableSetupColumnEx ("",
                            ImGuiTableColumnFlags_WidthFixed,
                            30.f, 0);

      IgTableSetupColumnEx ("Draw",
                            ImGuiTableColumnFlags_WidthFixed,
                            100.f, 0);

      IgTableSetupColumnEx ("Summary",
                            ImGuiTableColumnFlags_WidthStretch,
                            0.f, 0);

      for (int i = 0; i < draw_list.nitems; i++) {
        UIDrawItem *item = draw_list.items[i];
        bool selected = item == selected_draw_item;

        IgPushIDInt (i);

        IgTableNextRow ();

        if (selected)
          IgTableSetBgColor (ImGuiTableBgTarget_RowBg1,
                             IM_COL32 (80, 80, 120, 255), -1);
        
        /* Visibility */
        IgTableNextColumn ();
        IgBeginDisabled(item == view_item);
        if (IgCheckbox ("##visible", &item->visible)) {
          display_item_visible (item->id, item->visible);
          changed = true;
          if (item->visible)
            selected_draw_item = item;
          else if (selected_draw_item == item)
            selected_draw_item = NULL;
        }
        IgEndDisabled();
        
        /* Name */
        IgTableNextColumn ();
        IgPushStyleColor (ImGuiCol_HeaderHovered, IM_COL32 (0, 0, 0, 0));
        IgPushStyleColor (ImGuiCol_HeaderActive, IM_COL32 (0, 0, 0, 0));
        if (item->error)
          IgPushStyleColor (ImGuiCol_Text, IM_COL32(255, 0, 0, 255));
        if (IgSelectableEx (item->name, false, 0, (ImVec2){0, 0})) {
          selected_draw_item = item;
          properties_open = true;
        }
        IgPopStyleColorEx (item->error ? 3 : 2);
 
        /* Summary */
        IgTableNextColumn ();
        IgPushStyleColor (ImGuiCol_HeaderHovered, IM_COL32 (0, 0, 0, 0));
        IgPushStyleColor (ImGuiCol_HeaderActive, IM_COL32 (0, 0, 0, 0));
        if (item->error)
          IgPushStyleColor (ImGuiCol_Text, IM_COL32(255, 0, 0, 255));
        const char * summary = item->summary (item);
        if (IgSelectableEx (summary ? summary : "", false, 0, (ImVec2){0, 0})) {
          selected_draw_item = item;
          properties_open = true;
        }
        IgPopStyleColorEx (item->error ? 3 : 2);
 
        IgPopID ();
      }

      IgEndTable ();
    } // Draw items table

    /**
    ## Properties */
  
    if (ui_draw_properties ())
      changed = true;

    /**
    ## Status */

    if (UIState.status) {
      float y = IgGetContentRegionAvail().y;
      if (y > IgGetTextLineHeight())
        IgSetCursorPosY (IgGetCursorPosY() +
                         y - IgGetTextLineHeight());
      IgTextUnformatted (UIState.status);
    }
    
  } // Draw window

  IgEnd ();
  
  /**
  # Draggable separator. */

  const int separator_width = 10;
  
  IgPushStyleVarImVec2 (ImGuiStyleVar_WindowMinSize, (ImVec2){0.f, 0.f});
  IgSetNextWindowPos ((ImVec2){lwidth - separator_width/2.f, ui_top},
                      ImGuiCond_Always);
  IgSetNextWindowSize ((ImVec2){separator_width, height - ui_top},
                       ImGuiCond_Always);
  ImVec2 old_padding = IgGetStyle()->WindowPadding;
  IgGetStyle()->WindowPadding = (ImVec2){0.f, 0.f};
  if (IgBegin ("##draw_separator", NULL,
               ImGuiWindowFlags_NoDecoration |
               ImGuiWindowFlags_NoMove |
               ImGuiWindowFlags_NoResize |
               ImGuiWindowFlags_NoTitleBar |
               ImGuiWindowFlags_NoScrollbar |
               ImGuiWindowFlags_NoCollapse | 
               // ImGuiWindowFlags_NoSavedSettings |
               ImGuiWindowFlags_NoBringToFrontOnFocus)) {

    IgButtonEx ("##separator", IgGetContentRegionAvail());

    if (IgIsItemActivated ())
      draw_window_width_start = UIState.draw_window_width;

    if (IgIsItemActive ()) {
      ImVec2 delta = IgGetMouseDragDelta (0, -1.f);

      UIState.draw_window_width = draw_window_width_start + delta.x;

      if (UIState.draw_window_width < 150)
        UIState.draw_window_width = 150;
      if (UIState.draw_window_width > 600)
        UIState.draw_window_width = 600;
    }

  }

  IgEnd();
  
  IgPopStyleVar();
  IgGetStyle()->WindowPadding = old_padding;

  return changed;
}

/**
# Loading from file */

static UIDrawItem * process_line (const char * line)
{
  if (line[0] == '\0')
    return NULL;
  
  char copy[strlen(line) + 1];
  strcpy (copy, line);
  char * s = mystrtok (remove_blanks (copy), "(");
  if (!s)
    return NULL;

  if (!strcmp (s, "load")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      display_load (file);
    return NULL;
  }

  /**
  The [/src/draw_get.h]() file is generated automatically by [/src/params.awk]()
  and contains parsing commands for the functions defined in
  [/src/draw.h](). */

  #include "../../draw_get.h"

  return NULL;
}

bool ui_add_item (const char * command)
{
  UIDrawItem * item = process_line (command);
  if (item) {
    if (!strcmp (item->name, "view"))
      free (item);
    else
      ui_add_draw_item (item);
    return true;
  }
  return false;
}

/**
# Draw the user interface.

The UI belongs exclusively to the display thread. */

#include "bindings/dcimgui.h"
#include <GLFW/glfw3.h>

bool ui_display (GLFWwindow *window)
{
  bool changed = false;
  
  cImGui_ImplOpenGL3_NewFrame();
  cImGui_ImplGlfw_NewFrame();
  IgNewFrame();

  int width, height;
  
  glfwGetWindowSize (window, &width, &height);

  if (IgBeginMainMenuBar ()) {
    if (IgBeginMenu ("File")) {
#ifdef __EMSCRIPTEN__
      if (IgMenuItem ("Save"))
        gview_save ("display.gv");

      if (IgMenuItem ("Load"))
        EM_ASM ({ gview_open_file (); });
#else // !__EMSCRIPTEN__
      const char *filter[] = { "*.gv" };

      if (IgMenuItem ("Save")) {
        const char *filename = tinyfd_saveFileDialog ("Save gview", "display.gv",
                                                      1, filter, "gview files");
        if (filename)
          gview_save (filename);
      }
      
      if (IgMenuItem ("Load")) {
        const char *filename = tinyfd_openFileDialog ("Load gview", "",
                                                      1, filter, "gview files", 0);
        if (filename && display_load (filename))
          changed = true;
      }

      if (IgMenuItem ("Quit"))
        glfwSetWindowShouldClose (window, GLFW_TRUE);
#endif // !__EMSCRIPTEN__

      IgEndMenu ();
    }
    
    if (IgBeginMenu ("Edit")) {

      if (IgMenuItemEx ("Delete", NULL, false,
                        selected_draw_item != NULL &&
                        selected_draw_item != view_item)) {
        ui_delete_draw_item (selected_draw_item);
        changed = true;
      }
      
      if (IgMenuItemEx ("Clear", NULL, false, draw_list.nitems > 1)) {
        while (draw_list.nitems > 1) // fixme: always false!!!
          ui_delete_draw_item (*(draw_list.items + 1));
        changed = true;
      }

      IgEndMenu ();
    }
    
    if (IgBeginMenu ("View")) {
      if (UIState.dimension > 2) {
        if (IgMenuItemEx ("Front", "F", false, true))
          display_front_view();
        if (IgMenuItemEx ("Back",  "B", false, true))
          display_back_view();
        if (IgMenuItemEx ("Left", "L", false, true))
          display_left_view();
        if (IgMenuItemEx ("Right", "R", false, true))
          display_right_view();
        if (IgMenuItemEx ("Top", "T", false, true))
          display_top_view();
        if (IgMenuItemEx ("Bottom", "Ctrl+B", false, true))
          display_bottom_view();
        if (IgMenuItemEx ("Iso", "I", false, true))
          display_iso_view();
      }
      if (IgMenuItemEx ("Reset", "Ctrl+R", false, true))
        display_reset_view();
      IgEndMenu ();
    }

    if (ui_draw_menu())
      changed = true;

    IgEndMainMenuBar ();
  }

  if (ui_draw_list (width, height))
    changed = true;
  
  IgRender ();
  cImGui_ImplOpenGL3_RenderDrawData (IgGetDrawData ());

  return changed;
}

bool ui_want_capture_mouse()
{
  return IgGetIO()->WantCaptureMouse;
}

bool ui_want_capture_keyboard()
{
  return IgGetIO()->WantCaptureKeyboard;
}

void ui_viewport (GLFWwindow *window, int * xpos, int * ypos, int * width, int * height)
{
  glfwGetWindowSize (window, width, height);
  *xpos = UIState.draw_window_width;
  *ypos = 0;
  *width -= UIState.draw_window_width;
  *height -= ui_top;
  ui_view * ui = (ui_view *) view_item;
  ui->width = *width;
  ui->height = *height;
}

void ui_set_window_size (GLFWwindow *window, int width, int height)
{
  glfwSetWindowSize (window, width + UIState.draw_window_width, height + ui_top);
  ui_view *ui = (ui_view *) view_item;
  ui->width = width;
  ui->height = height;
}

/**
# State and display controls */

void ui_set_running (int running)
{
  UIState.running = running;
}

void ui_set_connected (int connected)
{
  if (UIState.connected == 1 && connected != 1)
    destroy_controls();
  UIState.connected = connected;
}

/**
# Initialisation and termination */

void ui_init (GLFWwindow *window)
{
  UIState.window = window;
  IgCreateContext (NULL);
  ImGuiIO *io = IgGetIO();
  ImFontAtlas_ClearFonts (io->Fonts);
  ImFontAtlas_AddFontFromFileTTF (io->Fonts,
#ifdef __EMSCRIPTEN__
                                  "/fonts/DejaVuSans.ttf",
#else
                                  "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
#endif // !__EMSCRIPTEN__
                                  font_size, NULL, NULL);
  cImGui_ImplGlfw_InitForOpenGL (window, true);
  cImGui_ImplOpenGL3_Init();

  /**
  The window size "memory" lives in `$HOME/.gview.ini`. */
     
  static char ini_filename[1024];
  snprintf (ini_filename, sizeof ini_filename,
            "%s/.gview.ini", getenv ("HOME"));
  IgGetIO ()->IniFilename = ini_filename;

  view_item = ui_view_new();
  ui_add_draw_item (view_item);
}

void ui_destroy()
{
  free (UIState.status);
  destroy_controls();
  cImGui_ImplOpenGL3_Shutdown();
  cImGui_ImplGlfw_Shutdown();
  IgDestroyContext (NULL);
}

