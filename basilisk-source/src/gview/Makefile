CFLAGS += -I$(BASILISK) -O2

# these are not Basilisk programs
EXCLUDE = display.c transfer_ws.c transfer_local.c wsclient.c

TOPTARGETS = all clean

SUBDIRS = gui

.PHONY: subdirs $(SUBDIRS) $(TOPTARGETS)

all: subdirs gview-nogui gview libdisplay.a transfer_local.o

subdirs: $(SUBDIRS)

$(TOPTARGETS): $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@ $(filter-out $(SUBDIRS),$(MAKECMDGOALS))

HEADERS = glad.h ../gl/trackball.h ../parse.h transfer.h display.h

display-nogui.o: display.c $(HEADERS)
	$(CC) -c $(CFLAGS) -o $@ $<

display.o: display.c $(HEADERS) gui/ui.h
	$(CC) -c $(CFLAGS) -DUI -o $@ $<

libdisplay.a: display.o glad.o \
	gui/ui.o gui/dcimgui_impl_glfw.o gui/dcimgui_impl_opengl3.o gui/dcimgui.o \
	gui/tinyfiledialogs.o
	ar cr $@ $^

gview: gview.o transfer_ws.o libdisplay.a
	$(CC) gview.o transfer_ws.o -o $@ -L. -ldisplay \
		-L../gl -lglutils -lglfw -limgui -lstb -lstdc++ -lm

gview-nogui: gview.o display-nogui.o transfer_ws.o glad.o
	$(CC) $^ -o $@ -L../gl -lglutils -lglfw -lm

wsclient: wsclient.o transfer_ws.o
	$(CC) $^ -o $@

wsclient.o: transfer.h

transfer_ws.o: transfer.h

transfer_local.o: transfer.h

EMCC_FLAGS = \
	-s USE_GLFW=3 \
	-s FULL_ES2=1 \
	-s ALLOW_MEMORY_GROWTH=1 \
	-s MAX_WEBGL_VERSION=2 \
	-sMINIFY_HTML=0 \
	-sEXIT_RUNTIME=1 \
	-sEXPORTED_RUNTIME_METHODS=ccall \
	--preload-file /usr/share/fonts/truetype/dejavu/DejaVuSans.ttf@/fonts/DejaVuSans.ttf

SRCS = gview.c display.c transfer_emscripten.c glad.c ../gl/trackball.c \
	gui/ui.c \
	gui/bindings/dcimgui_impl_glfw.cpp \
	gui/bindings/dcimgui_impl_opengl3.cpp \
        gui/bindings/dcimgui.cpp \
        gui/bindings/imgui.cpp \
        gui/bindings/imgui_widgets.cpp \
        gui/bindings/imgui_tables.cpp \
        gui/bindings/imgui_draw.cpp \
        gui/bindings/imgui_impl_glfw.cpp \
        gui/bindings/imgui_impl_opengl3.cpp

gview.js: $(SRCS) gui/ui.h
	emcc -Wall -I$(BASILISK) -I/usr/include/imgui -I/usr/include/imgui/backends -O2 $(EMCC_FLAGS) -DUI -Wno-unused-function -Wno-format-security $(SRCS) -o $@ -lwebsocket -lm

clean_emcc:
	rm -f gview.wasm gview.js gview.data gview.wasm.map

clean: clean_emcc

include ../Makefile.defs
