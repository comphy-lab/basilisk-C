TEMPLATES=$(BASILISK)/darcsit/templates/

all: page.static edit.static

%.static: *.st $(TEMPLATES)includes.awk
	awk -v path=$(TEMPLATES) -f $(TEMPLATES)includes.awk $*.st | \
	awk -v path=$(TEMPLATES) -f $(TEMPLATES)includes.awk | \
	awk -v path=$(TEMPLATES) -f $(TEMPLATES)includes.awk | \
	sed 's/\\\$$/__ESCAPEDDOLLAR__/g' > $@
