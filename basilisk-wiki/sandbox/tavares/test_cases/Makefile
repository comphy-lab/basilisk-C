include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in ../config
CFLAGS += -O2 -DMTRACE=3

impact-thin-fiber2.c: impact-thin-fiber.c
	ln -s impact-thin-fiber.c impact-thin-fiber2.c
impact-thin-fiber2.s: CFLAGS += -DCASE2=1
impact-thin-fiber2.tst: CFLAGS += -DCASE2=1

impact-thin-fiber.tst: impact-thin-fiber2.tst