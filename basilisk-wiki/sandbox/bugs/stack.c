/**
# Looks like a stack overflow */

#include "view.h"

int main() {
  for (int i = 0; i < 90; i++) {
    view (camera = "front", width = 800, height = 800);
    save ("a.ppm");
  }
}
