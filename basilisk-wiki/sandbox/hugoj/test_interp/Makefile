

all: test_interp test_interp_pol2cart


interpolate.c:
	ln -sf ${MYSANDBOX}/hugoj/lib/interpolate.h interpolate.h

test_interp: test_interp.c interpolate.c
	#gcc -o test_interp test_interp.c interpolate.c -lm
	gcc -o test_interp test_interp.c -lm
	./test_interp

test_interp_pol2cart: test_interp_pol2cart.c interpolate.c
	gcc -o test_interp_pol2cart test_interp_pol2cart.c -lm
	./test_interp_pol2cart


clean:
	rm -f test_interp test_interp_pol2cart

