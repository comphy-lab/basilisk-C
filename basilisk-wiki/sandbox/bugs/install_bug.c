/**
# Basilisk installation bug
If the username of the user starts with ast, then an error occurs during the installation steps of basilisk C. Namely, in the -DBASILISK compilation flag, the /ast is skipped, shortening e.g.: /home/astroman/basilisk/src\" to /homeroman/basilisk/src\". This behaviour is not seen if the below steps are made for the user fastroman. 

## Steps to reproduce
1. First create a new user starting with ast, using e.g.
   - sudo useradd --create-home astroman
   - sudo passwd astroman
2. Log into new user:
   - su - astroman
3. Install basilisk:
   - wget https://basilisk.fr/basilisk/basilisk.tar.gz
   - tar xzf basilisk.tar.gz
   - cd basilisk/src/
   - ln -s config.gcc config
   - make -k

## Example output from bug

sh /home/astroman/basilisk/src/tests.sh

updating /home/astroman/basilisk/src/Makefile.tests

updating Makefile.deps

make -C darcsit 

make[1]: Entering directory '/home/astroman/basilisk/src/darcsit'

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DYY_NO_UNPUT    literate-c.c   -o literate-c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DYY_NO_UNPUT    codeblock.c   -o codeblock

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DYY_NO_UNPUT    pagemagic.c   -o pagemagic

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DYY_NO_UNPUT -DYY_NO_INPUT    sanitize.c   -o sanitize

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DYY_NO_UNPUT    urldecode.c   -o urldecode

make[1]: Leaving directory '/home/astroman/basilisk/src/darcsit'

make -C ast 

make[1]: Entering directory '/home/astroman/basilisk/src/ast'

make -C interpreter 

make[2]: Entering directory '/home/astroman/basilisk/src/ast/interpreter'

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\" -I../.. -I.. -DBASILISK=\"/home/astroman/basilisk/src\"   -c -o interpreter.o interpreter.c

<command-line>: warning: ‘BASILISK’ redefined

<command-line>: note: this is the location of the previous definition

make[2]: Leaving directory '/home/astroman/basilisk/src/ast/interpreter'

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\"   -c -o ast.o ast.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\"   -c -o tokens.o tokens.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\"   -c -o basilisk.o basilisk.c

basilisk.c: In function ‘yyparse’:

basilisk.c:4031:9: warning: variable ‘yynerrs’ set but not used [-Wunused-but-set-variable=]

 4031 |     int yynerrs = 0;
 
      |         ^~~~~~~
      
cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\"   -c -o translate.o translate.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\"   -c -o allocator.o allocator.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\"   -c -o faststack.o faststack.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\"   -c -o stencil.o stencil.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\"   -c -o types.o types.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\"   -c -o references.o references.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\"   -c -o kernels.o kernels.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\" grammar.c -o grammar

./grammar < basilisk.yacc  > grammar.h

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -DBASILISK=\"/homeroman/basilisk/src\"   -c -o check.o check.c

ar cr libast.a ast.o tokens.o basilisk.o translate.o allocator.o faststack.o

stencil.o types.o references.o kernels.o check.o interpreter/interpreter.o

make[1]: Leaving directory '/home/astroman/basilisk/src/ast'

make -C kdt 

make[1]: Entering directory '/home/astroman/basilisk/src/kdt'

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -D_FILE_OFFSET_BITS=64 -c kdt.c

ar cr libkdt.a kdt.o

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 

xyz2kdt.c kdt.o -o xyz2kdt -lm

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2
kdtquery.c kdt.o -o kdtquery -lm

make[1]: Leaving directory '/home/astroman/basilisk/src/kdt'

make -C wsServer 

make[1]: Entering directory '/home/astroman/basilisk/src/wsServer'

cc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1 -pedantic  -Wno-unused-result -Wno-overlength-strings -fno-diagnostics-show-caret /home/astroman/basilisk/src/wsServer/src/base64/base64.c -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -I /home/astroman/basilisk/src/wsServer/include -c -o /home/astroman/basilisk/src/wsServer/src/base64/base64.o

cc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1 -pedantic  -Wno-unused-result -Wno-overlength-strings -fno-diagnostics-show-caret /home/astroman/basilisk/src/wsServer/src/handshake/handshake.c -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -I /home/astroman/basilisk/src/wsServer/include -c -o /home/astroman/basilisk/src/wsServer/src/handshake/handshake.o

cc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1 -pedantic  -Wno-unused-result -Wno-overlength-strings -fno-diagnostics-show-caret /home/astroman/basilisk/src/wsServer/src/sha1/sha1.c -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -I /home/astroman/basilisk/src/wsServer/include -c -o /home/astroman/basilisk/src/wsServer/src/sha1/sha1.o

cc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1 -pedantic  -Wno-unused-result -Wno-overlength-strings -fno-diagnostics-show-caret /home/astroman/basilisk/src/wsServer/src/ws.c -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -I /home/astroman/basilisk/src/wsServer/include -c -o /home/astroman/basilisk/src/wsServer/src/ws.o

ar cr libws.a /home/astroman/basilisk/src/wsServer/src/base64/base64.o /home/astroman/basilisk/src/wsServer/src/handshake/handshake.o /home/astroman/basilisk/src/wsServer/src/sha1/sha1.o /home/astroman/basilisk/src/wsServer/src/ws.o

make[1]: Leaving directory '/home/astroman/basilisk/src/wsServer'

make -C gl 

make[1]: Entering directory '/home/astroman/basilisk/src/gl'

make -C tinyrenderer 

make[2]: Entering directory '/home/astroman/basilisk/src/gl/tinyrenderer'

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -c tiny.c

make[2]: Leaving directory '/home/astroman/basilisk/src/gl/tinyrenderer'

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I..   -c -o trackball.o trackball.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I..   -c -o utils.o utils.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I..   -c -o polygonize.o polygonize.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I..   -c -o og_stroke_mono_roman.o og_stroke_mono_roman.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I.. -c -o parser.o parser.tab.c

parser.tab.c: In function ‘yyparse’:

parser.tab.c:1041:9: warning: variable ‘yynerrs’ set but not used [-Wunused-but-set-variable=]

 1041 |     int yynerrs = 0;
 
      |         ^~~~~~~
      
cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I..   -c -o TinyPngOut.o TinyPngOut.c

TinyPngOut.c: In function ‘TinyPngOut_init’:

TinyPngOut.c:37:30: warning: ‘nonnull’ argument ‘out’ compared to NULL [-Wnonnull-compare]

   37 |         if (w == 0 || h == 0 || out == NULL)
     
      |                       
      ^
ar cr libglutils.a trackball.o utils.o polygonize.o og_stroke_mono_roman.o parser.o TinyPngOut.o

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -I..   -c -o fb_tiny.o fb_tiny.c

ar cr libfb_tiny.a fb_tiny.o tinyrenderer/tiny.o

make[1]: Leaving directory '/home/astroman/basilisk/src/gl'

make -C grid 

make[1]: Entering directory '/home/astroman/basilisk/src/grid'

make -C gpu 

make[2]: Entering directory '/home/astroman/basilisk/src/grid/gpu'

sh /home/astroman/basilisk/src/tests.sh

updating /home/astroman/basilisk/src/grid/gpu/Makefile.tests

updating Makefile.deps

make[2]: *** No rule to make target 'errors.c.tags.d', needed by 'Makefile.deps'.

make[2]: *** No rule to make target 'glad.c.tags.d', needed by 'Makefile.deps'.

make[2]: *** No rule to make target 'glad.h.tags.d', needed by 'Makefile.deps'.

make[2]: *** No rule to make target 'grid.h.tags.d', needed by 'Makefile.deps'.

make[2]: *** No rule to make target 'opengl.c.tags.d', needed by 'Makefile.deps'.

make[2]: *** No rule to make target 'output.h.tags.d', needed by 'Makefile.deps'.

../../Makefile.defs:248: Failed to remake makefile 'Makefile.deps'.

make[2]: Leaving directory '/home/astroman/basilisk/src/grid/gpu'

make[1]: *** [Makefile:11: gpu] Error 2
  
make -C cuda 

make[2]: Entering directory '/home/astroman/basilisk/src/grid/cuda'

sh /home/astroman/basilisk/src/tests.sh

updating /home/astroman/basilisk/src/grid/cuda/Makefile.tests

updating Makefile.deps

make[2]: *** No rule to make target 'cuda.c.tags.d', needed by 'Makefile.deps'.
../../Makefile.defs:248: Failed to remake makefile 'Makefile.deps'.

make[2]: Leaving directory '/home/astroman/basilisk/src/grid/cuda'

make[1]: *** [Makefile:11: cuda] Error 2

make -C hip 

make[2]: Entering directory '/home/astroman/basilisk/src/grid/hip'
  
sh /home/astroman/basilisk/src/tests.sh

updating /home/astroman/basilisk/src/grid/hip/Makefile.tests

updating Makefile.deps

make[2]: *** No rule to make target 'cartesian.h.tags.d', needed by 'Makefile.deps'.

make[2]: *** No rule to make target 'hip.c.tags.d', needed by 'Makefile.deps'.

make[2]: *** No rule to make target 'multigrid.h.tags.d', needed by 'Makefile.deps'.

../../Makefile.defs:248: Failed to remake makefile 'Makefile.deps'.

make[2]: Leaving directory '/home/astroman/basilisk/src/grid/hip'

make[1]: *** [Makefile:11: hip] Error 2
  
make -C opencl 

make[2]: Entering directory '/home/astroman/basilisk/src/grid/opencl'

sh /home/astroman/basilisk/src/tests.sh

updating /home/astroman/basilisk/src/grid/opencl/Makefile.tests

updating Makefile.deps

make[2]: *** No rule to make target 'opencl.c.tags.d', needed by 'Makefile.deps'.

../../Makefile.defs:248: Failed to remake makefile 'Makefile.deps'.

make[2]: Leaving directory '/home/astroman/basilisk/src/grid/opencl'

make[1]: *** [Makefile:11: opencl] Error 2
  
make[1]: Target 'subdirs' not remade because of errors.

make[1]: Leaving directory '/home/astroman/basilisk/src/grid'

make: *** [Makefile:23: grid] Error 2

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DLIBDIR=\"`pwd`\" -c include.c
  
cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DLIBDIR=\"`pwd`\" -c postproc.c

cc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DLIBDIR=\"`pwd`\" \
        -DCC99="\"cc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1 -pedantic  -Wno-unused-result -Wno-overlength-strings -fno-diagnostics-show-caret\"" \
        -DCPP99="\"\"" \
        -DCADNACC="\"clang -D_CADNA=1 -x c++ -m64 -Wno-unused-function -Wno-unused-result -Wno-c++11-compat-deprecated-writable-strings -Wno-address-of-array-temporary\"" \
        -DBASILISK="\"/home/astroman/basilisk/src\"" \
        qcc.c include.o postproc.o -o qcc -Last -last -lm
 
cd darcsit && make && cd cgi-bin && make

make[1]: Entering directory '/home/astroman/basilisk/src/darcsit'

make[1]: Nothing to be done for 'all'.

make[1]: Leaving directory '/home/astroman/basilisk/src/darcsit'

make[1]: Entering directory '/home/astroman/basilisk/src/darcsit/cgi-bin'

chmod 755 darcsit

make[1]: Leaving directory '/home/astroman/basilisk/src/darcsit/cgi-bin'

./qcc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -autolink bview.c -o bview2D -lfb_tiny -lm

qcc: translate.c:5061: endfor: Assertion `fp' failed.

make: *** [Makefile:58: bview2D] Avbrutt (SIGABRT) (core dumped)

./qcc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -autolink -grid=octree bview.c -o bview3D -lfb_tiny -lm

qcc: translate.c:5061: endfor: Assertion `fp' failed.

make: *** [Makefile:61: bview3D] Avbrutt (SIGABRT) (core dumped)

./qcc -std=c99 -D_XOPEN_SOURCE=700 -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -autolink -grid=multigrid bview.c -o bview2Dm -lfb_tiny -lm

qcc: translate.c:5061: endfor: Assertion `fp' failed.

make: *** [Makefile:64: bview2Dm] Avbrutt (SIGABRT) (core dumped)

make: Target 'all' not remade because of errors.

*/