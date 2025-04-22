# if you want to use a memory-buffered KDT, add -DBUFFERED_KDT=1 to CFLAGS
CFLAGS += -O2

all: libkdt.a xyz2kdt kdtquery

libkdt.a: kdt.o
	ar cr libkdt.a kdt.o

kdt.o: kdt.c
	$(CC) $(CFLAGS) -D_FILE_OFFSET_BITS=64 -c kdt.c

xyz2kdt: xyz2kdt.c kdt.o kdt.h
	$(CC) $(CFLAGS) xyz2kdt.c kdt.o -o xyz2kdt -lm

kdtquery: kdtquery.c kdt.o kdt.h
	$(CC) $(CFLAGS) kdtquery.c kdt.o -o kdtquery -lm

clean:
	rm -f *.o

check:
