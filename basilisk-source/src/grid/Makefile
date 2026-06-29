TOPTARGETS = all clean check

SUBDIRS = gpu cuda hip opencl

.PHONY: subdirs $(SUBDIRS) $(TOPTARGETS)

subdirs: $(SUBDIRS)

$(TOPTARGETS): $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@ $(filter-out $(SUBDIRS),$(MAKECMDGOALS))

all:
clean:
check:
