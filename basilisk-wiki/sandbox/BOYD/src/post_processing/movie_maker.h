/**
# Movie maker that adds scalar limits
*/
double TX = -0.5;
double TY = -0.5;
double TZ = 0.0;
double FOV = 24;
float QUAT[4] = {0, 0, 0, 0};
int MESH_ON = 1;
int BOX_ON = 0;
int SHOW_LIMITS = 1;
int SHOW_TIME = 1;
bool LINEAR_ = true;
double THETA_CAMERA = 0.0;
float NORM[3] = {0, 0, 1};
double ALPHA_ = 0.0;
bool RELATIVE_ = false;
bool MUTE_VIEW = false;
float NEAR_ = 0.0;
float FAR_ = 0.0;

#ifndef BG
#define BG 0.7  // light gray for background
#endif
#ifndef DG
#define DG 0.   // dark gray
#endif

double scalar_max(scalar P) {
  double P_max = -1.0e12;
  foreach (reduction(max:P_max)){
    if (P[] > P_max) {
      P_max = P[];
    }
  }
  return P_max;
}
//
double scalar_min(scalar P) {
  double P_min = 1.0e12;
  foreach (reduction(min:P_min)){
    if (P[] < P_min) {
      P_min = P[];
    }
  }
  return P_min;
}

void movie_maker_i(char* var_name, int i, double t) {
  scalar P = lookup_field(var_name);
  double P_min = scalar_min(P);
  double P_max = scalar_max(P);
  //
  if (!MUTE_VIEW) {
    view(fov = FOV, quat = {QUAT[0], QUAT[1], QUAT[2], QUAT[3]},
       width = 640, height = 640, relative = RELATIVE_,
       tx = TX, ty = TY, bg = {BG, BG, BG}, view = NULL);
    MUTE_VIEW = false;
  }
  if (MESH_ON) {
    cells();
  }
  if (BOX_ON) {
    box();
  }

  squares(var_name, min = P_min, max = P_max, linear = LINEAR_,
          map = cool_warm);

  draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0);
  char str[99];
  sprintf(str, "  %10g < %s < %10g", P_min, var_name, P_max);
  if (SHOW_LIMITS) {
    draw_string(str, size = 60, lw = 2, pos = 1);
  }
  char t_str[99];
  sprintf(t_str, "t = %10g", t);
  if (SHOW_TIME) {
    draw_string(t_str, size = 60, lw = 2, pos = 2);
  }

  const char* ext = ".mp4";

  strcpy(str, var_name);
  strcat(str, ext);
  save(str);
}
#if dimension == 3
#include "iso3D.h"
void vof_movie_3D_mirror_x_y_z() {
  if (!MUTE_VIEW) {
    view(camera = "iso", fov = FOV, width = 640, height = 640,
         relative = RELATIVE_, tx = TX, ty = TY, bg = {BG, BG, BG});
    MUTE_VIEW = false;
  }
  draw_vof("f");
  if (MESH_ON) {
    cells();
  }
  if (BOX_ON) {
    box();
  }
  mirror(n = {1, 0, 0}) {
    draw_vof("f");
    if (MESH_ON) {
      cells();
    }
  }
  mirror(n = {0, 1, 0}) {
    draw_vof("f");
    if (MESH_ON) {
      cells();
    }
    mirror(n = {1, 0, 0}) {
      draw_vof("f");
      if (MESH_ON) {
        cells();
      }
    }
  }
  mirror(n = {0, 0, 1}) {
    draw_vof("f");
    if (MESH_ON) {
      cells();
    }
    mirror(n = {1, 0, 0}) {
      draw_vof("f");
      if (MESH_ON) {
        cells();
      }
    }
    mirror(n = {0, 1, 0}) {
      draw_vof("f");
      if (MESH_ON) {
        cells();
      }
      mirror(n = {1, 0, 0}) {
        draw_vof("f");
        if (MESH_ON) {
          cells();
        }
      }
    }
  }
  box();
  save("vof_movie.mp4");
}

void movie_maker_3D_mirror_x_y(char* var_name, char* pov, double t) {
  scalar P = lookup_field(var_name);
  double P_min = scalar_min(P);
  double P_max = scalar_max(P);
  if (!MUTE_VIEW) {
    view(camera = pov, fov = FOV, width = 640, height = 640,
         relative = RELATIVE_, tx = TX, ty = TY, bg = {BG, BG, BG});
    MUTE_VIEW = false;
  }

  cross_section("f", alpha = 0.0);
  if (MESH_ON) {
    cells();
  }
  if (BOX_ON) {
    box();
  }
  squares(var_name, n = {NORM[0], NORM[1], NORM[2]}, alpha = ALPHA_,
        min = P_min, max = P_max, linear = LINEAR_, map = cool_warm);
  mirror(n = {1, 0, 0}) {
    cross_section("f", alpha = 0.0);
    if (MESH_ON) {
      cells();
    }
    squares(var_name, n = {NORM[0], NORM[1], NORM[2]}, alpha = ALPHA_,
          min = P_min, max = P_max, linear = LINEAR_, map = cool_warm);
  }
  mirror(n = {0, 1, 0}) {
    cross_section("f", alpha = 0.0);
    if (MESH_ON) {
      cells();
    }
    squares(var_name, n = {NORM[0], NORM[1], NORM[2]}, alpha = ALPHA_,
            min = P_min, max = P_max, linear = LINEAR_, map = cool_warm);
    mirror(n = {1, 0, 0}) {
      cross_section("f", alpha = 0.0);
      if (MESH_ON) {
        cells();
      }
      squares(var_name, n = {NORM[0], NORM[1], NORM[2]}, alpha = ALPHA_,
              min = P_min, max = P_max, linear = LINEAR_,
              map = cool_warm);
    }
  }

  char str[99];
  sprintf(str, "  %10g < %s < %10g", P_min, var_name, P_max);
  if (SHOW_LIMITS) {
    draw_string(str, size = 60, lw = 2, pos = 1);
  }
  char t_str[99];
  sprintf(t_str, "t = %10g", t);
  if (SHOW_TIME) {
    draw_string(t_str, size = 60, lw = 2, pos = 2);
  }

  const char* ext = ".mp4";

  strcpy(str, var_name);
  strcat(str, ext);
  save(str);
}
#endif

