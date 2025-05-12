#include "view.h"

int main(){
  init_grid (N);
        view();
        box();
        save("box.png");
}

/**
> is this a bug or what did i do wrong?

There needs to be a grid...

![Img](test_box/box.png)


1339 /home/basilisk-src/basilisk/src/draw.h: No such file or directory.
src/draw.h:1339:error: Program received signal SIGSEGV, Segmentation fault.
*/