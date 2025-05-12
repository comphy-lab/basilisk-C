/**
# Test the puzzle generator

![a hard puzzle?](test_puzzle/puzzle.mp4)
*/

#include "view.h"
#include "puzzle.h"

int main() {
  view (tx = -0.5, ty = -0.5);
  for (int frame = 20; frame < 100; frame++) {
    init_grid (4);
    refine (x + y < frame/50. && level < 6);
    unrefine (x + y < frame/50. - 0.1);
    puzzle_cells(lw = 2);
    save ("puzzle.mp4");
  }
}
