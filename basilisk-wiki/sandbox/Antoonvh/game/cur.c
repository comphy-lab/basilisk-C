#include <ncurses.h>

void start_screen(void) {
  initscr();
  noecho();
  cbreak();
  nodelay(stdscr, true);
}

void stop_screen(void) {
  endwin();
}
