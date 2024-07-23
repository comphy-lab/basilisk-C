CFLAGS += $(C99_FLAGS) -O2

export BASILISK = $(CURDIR)

# these are not Basilisk programs
EXCLUDE = qcc.c include.c postproc.c bview.c

TOPTARGETS = all clean check

SUBDIRS = darcsit ast kdt wsServer gl

.PHONY: subdirs $(SUBDIRS) $(TOPTARGETS)

all: subdirs qcc literatec bview2D bview3D
	@chmod +x ppm2mpeg ppm2mp4 ppm2ogv ppm2gif runtest page2html
	@test -f xyz2kdt || ln -s kdt/xyz2kdt
	@test -f kdtquery || ln -s kdt/kdtquery

subdirs: $(SUBDIRS)

$(TOPTARGETS): $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@ $(filter-out $(SUBDIRS),$(MAKECMDGOALS))

literatec:
	cd darcsit && $(MAKE) && cd cgi-bin && $(MAKE)

qcc: qcc.c include.o postproc.o ast/libast.a config
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" \
		-DCC99="\"$(CC99)\"" \
		-DCPP99="\"$(CPP99)\"" \
		-DCADNACC="\"$(CADNACC)\"" \
		-DBASILISK="\"$(BASILISK)\"" \
		qcc.c include.o postproc.o -o qcc -Last -last -lm

include.o: include.c
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" -c include.c

postproc.o: postproc.c
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" -c postproc.c

# uncomment the recipe below if you need to regenerate draw_get.h
# and draw_json.h

# draw_get.h: draw.h params.awk
#	awk -f params.awk < draw.h > draw_get.h

# Uncomment the recipes below if you need to re-generate include.c or
# postproc.c

# include.c: include.lex
#	flex -P inc -o include.c include.lex

# postproc.c: postproc.lex
#	flex -P post -o postproc.c postproc.lex

bview2D: qcc bview.c draw_get.h display.h view.h draw.h khash.h
	./qcc $(CFLAGS) -autolink bview.c -o bview2D -lfb_tiny -lm

bview3D: qcc bview.c draw_get.h display.h view.h draw.h khash.h
	./qcc $(CFLAGS) -autolink -grid=octree bview.c -o bview3D -lfb_tiny -lm

alltags:
	cd navier-stokes && $(MAKE) tags
	cd layered && $(MAKE) tags
	cd ehd && $(MAKE) tags
	cd compressible && $(MAKE) tags
	cd examples && $(MAKE) tags
	cd test && $(MAKE) tags
	$(MAKE) tags
	cd navier-stokes && $(MAKE) itags
	cd layered && $(MAKE) itags
	cd compressible && $(MAKE) itags
	cd ehd && $(MAKE) itags
	cd examples && $(MAKE) itags
	cd test && $(MAKE) itags
	$(MAKE) itags

etags:
	etags *.h grid/*.h

checklinks:
	$(LINKCHECKER) 	$(BASILISK_URL)/src/README 		\
			$(BASILISK_URL)/src/test/README 	\
			$(BASILISK_URL)/src/examples/README | 	\
		tee checklinks.log

checklinksfast:
	wget --spider -nd -nv -r \
		--reject-regex '.*[?]changes=.*' \
		--reject-regex '.*[?]history' \
		$(BASILISK_URL) 2>&1 | \
		grep -v ^unlink: | tee checklinks.log

changelog:
	darcs changes > ChangeLog

dist:
	darcs dist

diff:
	cd .. && tar czvf src/diff.tgz `darcs whatsnew -s | \
		sed 's/. .\/.*\/$$//g' | awk '{print $$2}'`

include $(BASILISK)/Makefile.defs
