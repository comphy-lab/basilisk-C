# PYTHON_CONFIG = python3.10-config
# PYTHON_INCLUDE = $(shell $(PYTHON_CONFIG) --includes)
# PYTHON_LIBS = $(shell $(PYTHON_CONFIG) --ldflags)

#exec=test_interp
exec=test_interp_pol2cart

all: $(exec)
$(exec): $(exec).c interpolate.c
	#$(CC) -o $(exec) $(exec).c interpolate.c $(PYTHON_INCLUDE) $(PYTHON_LIBS)
	$(CC) -o $(exec) $(exec).c interpolate.c -lm

clean:
	rm -f $(exec)

