/**
# OSmesa package info on the sandbox server

see [here](OSMesainfo/out) and [here](OSMesainfo/log)
*/
#include "view.h"

int main() {
  system ("apt show -a libosmesa6-dev");
  system ("uname -a");
  system ("ls /home/basilisk/lib");
  system ("apt list --installed | grep osmesa 1>&2");
}