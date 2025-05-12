#include "common.h"

struct thatStruct {
  double a;
};

struct thatStruct var;

@define p() var

struct thisStruct {
  double x;
};

void func(struct thisStruct p) {
  p().a = 500;
}

int main() {
  var.a = 0;
  struct thisStruct p;
  p.x = 100;
  func(p);
  return 0;
}