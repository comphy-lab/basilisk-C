.PHONY: all clean

CFLAGS += -I/usr/include/imgui -O2

all: ui.o tinyfiledialogs.o dcimgui_impl_glfw.o dcimgui_impl_opengl3.o dcimgui.o

ui.o: ui.h ../display.h tinyfiledialogs.h ../../draw_ui.h ../../draw_get.h ../../parse.h

# these are not Basilisk programs
EXCLUDE = ui.c tinyfiledialogs.c

%.o: bindings/%.cpp
	g++ $(CFLAGS) -I/usr/include/imgui -I/usr/include/imgui/backends -c -o $@ $<

include ../../Makefile.defs
