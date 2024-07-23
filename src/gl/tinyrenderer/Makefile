tiny.o: tiny.c tiny.h geometry.h ../framebuffer.h
	$(CC) $(CFLAGS) -c tiny.c

# An example of usage, requires the models in
# https://github.com/ssloy/tinyrenderer/tree/master/obj
#
# make main
# ./main obj/diablo3_pose/diablo3_pose.obj

OBJS = model.cpp tgaimage.cpp main.cpp
HDS = geometry.h  model.h  tiny.h tgaimage.h

main: $(OBJS) tiny.o $(HDS)
	g++ $(CFLAGS) $(OBJS) tiny.o -o main

clean:
	rm -f tiny.o main
