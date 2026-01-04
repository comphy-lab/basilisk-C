# where is the GOTM source code?
ifeq ($(GOTM),)
	GOTM = $(HOME)/local/src/GOTM
endif

# matching directories in the GOTM source code $GOTM/src/
DIRS = airsea gotm input meanflow observations turbulence util

all:
	for d in $(DIRS); do ( cd $$d && make -f ../Makefile headers ) done

MODULE = $(notdir $(PWD))

HEADERS = $(subst $(GOTM)/src/$(MODULE)/,,	\
	$(subst .F90,.h,$(shell ls $(GOTM)/src/$(MODULE)/*.F90)))

headers: $(HEADERS)

%.h: $(GOTM)/src/$(MODULE)/%.F90 ../comments.awk ../gotm.awk
	@echo "extracting $@ from $<"
	@( tr '[:upper:]' '[:lower:]' < $< | awk -f ../comments.awk 	\
	| awk -f ../gotm.awk -vbase=$* -v sourcefile="$<" 		\
		-v dirname=$(MODULE) > $@ ) || 				\
	( rm -f $*.F90 && exit 1 )

dependencies.svg: comments.awk dependencies.awk
	cat `find $(GOTM)/src -name '*.F90'` | awk -f comments.awk | \
	awk -f dependencies.awk | unflatten | dot -Tsvg > dependencies.svg
