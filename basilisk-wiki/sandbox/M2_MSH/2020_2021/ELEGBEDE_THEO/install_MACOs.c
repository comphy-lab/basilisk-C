/**
# Installing Basilisk on Mac OS.

Installing basilisk itself is quiet trivial. The tricky parts come from the tools used to obtain visualisation of the results. We shall start with that.

### Bview installation

First, you'll need to check that you have Homebrew installed on your mac. Run the following command line :

*/

brew --version

/**
Which on my Mac returns :
*/

Homebrew 2.5.9
Homebrew/homebrew-core (git revision 2cd74; last commit 2020-11-10)

/**
If you get a result stating that the command was not found, then run the following command on your terminal :
*/

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

/**


Then, you'll need to install pkgconfig. Run :

*/
brew install pkgconfig

/**
To avoid having to write the whole path to the directory containging pkg, we'll create a variable that contains that path :
*/

export PKG_CONFIG_PATH="~/soft/mesa/lib/pkgconfig:/opt/X11/lib/pkgconfig"
  
/**
You will then need to install meson and ninja. Enter the two following lines : 
*/
  
brew install meson
brew install ninja
 
 /**
 
 You'll need to create a directory where you'll stock the directories and files you'll create later. Do :
 
*/

mkdir ~/soft

/**
 Once meson and ninja are installed, you'll need to havec access to the x11 servers. To do so, download XQuartz from the following link :
 
 Depuis le lien suivant, installer le serveur XQuartz sur macOS pour accéder a X11 : 
https://www.xquartz.org/

Once it is done, open a finder windows. Go to your user directory and press cmd + shift + dot to reveal the hidden files. Open the .bash_profile file, and at the end of that file, add :
 */
 
export PKG_CONFIG_PATH=/opt/X11/lib/pkgconfig
export PKG_CONFIG_PATH=~/soft/mesa/lib/pkgconfig:/opt/X11/lib/pkgconfig
 
 /**
 Now follow this link : https://archive.mesa3d.org// 
 At the bottom of the page, download the latest version of mesa. The name should be : mesa-XX.X.X.tar.xz, where the Xs stand for the version you'll choose.
 
 Remplacing mesa-XX.X.X.tar.xz by the name of your file, do : 
 */
 
mv ~/Downloads/mesa-XX.X.X.tar.xz ~/Documents/Work/mesa-XX.X.X.tar.xz
tar xzf mesa-XX.X.X.tar.xz

/**
Now a file mesa-XX.X.X should have appeared in your soft directory. Enter that directory and open the file meson.build. Modify the first section that starts with project so that it is identical to : 
*/

project(
  'mesa',
  ['c', 'cpp'],
  version : run_command(
    [find_program('python', 'python2', 'python3'), 'bin/meson_get_version.py']
  ).stdout(),
  license : 'MIT',
  meson_version : '>= 0.46',
  default_options : ['buildtype=debugoptimized', 'b_ndebug=if-release', 'c_std=c11', 'cpp_std=c++11']
)

/**
Now, you'll need a specific C++ compiler and to have a specific version of mako. Assuming that anaconda is already installed on your mac (if not follow that link : https://docs.anaconda.com/anaconda/install/mac-os/), in your terminal, do :
*/

brew install gcc@9
conda install mako

/**
Once it is installed :
*/

cd ~/soft/mesa-XX.X.X
CC="gcc-9" CXX="gcc-9" CFLAGS="-I/~/soft/mesa/include/ -I/opt/X11/include/ -I/usr/local/include/" CXXFLAGS="-I/~soft/mesa/include/ -I/opt/X11/include/ -I/usr/local/include/" LDFLAGS=" -L/opt/X11/lib/ -L/usr/local/lib" meson --prefix=~/soft/mesa/ -Dosmesa=gallium build/
 
/**
Then :
*/
  
meson configure build/ | grep osmesa
meson configure -Dosmesa=gallium -Dosmesa-bits=8  build/
  
/** Now comes a command that should run for quiet some time : 
*/
  
ninja -C build/
ninja -C build/ install

/**  In the .bash_profile : 
*/

export LD_LIBRARY_PATH=~/mesa/lib:$LD_LIBRARY_PATH

/** Then, you will need some more installations :
*/

brew install autoconf
brew install automake
brew install libtool

/* You will then need to access each of the 2 following websites to download the glu, glut and glew library : 
- https://gitlab.freedesktop.org/mesa/glu
- https://gitlab.freedesktop.org/mesa/glut

Once installed, go to each of the folders obtained, and do : */

CFLAGS="-I/~/soft/mesa/include/ -I/opt/X11/include/" LDFLAGS="-L/~/soft/mesa/lib/ -L/opt/X11/lib/" ./autogen.sh --prefix=~/soft/mesa

/** Once its done, in the glu directory do : */

make
make install

/** in the glut directory do : */ 

brew install makedepend
make
make install

/** Then for a third library : */


brew install glew
cd ~/basilisk/src
cp config.osx config

/** Then open the config text file you just created and at its end add : */

OPENGLIBS = -L/Users/Olaitan/soft/mesa/lib/ -lfb_osmesa -lGLU -lOSMesa

/** Back to your terminal : */ 

cd gl

/** Open the fb_osmesa.c text file and change the first line to : */

#include </Users/Olaitan/soft/mesa/include/GL/osmesa.h>

/** Once again, back to your terminal : */

CC="gcc-9" CFLAGS="-I/~/soft/mesa/include/ -I/opt/X11/include/ -I/usr/local/include/" make libglutils.a libfb_osmesa.a
brew install imagemagick
cd ..
  
/** Open the config text file and change, change it entirely to the following version : */

  # -*-Makefile-*-
CC=gcc-9
CXX=g++-9
CFLAGS += -I/~/soft/mesa/include/ -I/opt/X11/include/ -I/usr/local/include/
CXXFLAGS += -I/~soft/mesa/include/ -I/opt/X11/include/ -I/usr/local/include/
LDFLAGS += -L/~/soft/mesa/lib/ -L/opt/X11/lib/ -L/usr/local/lib

# how to launch the C99 compiler
CC99 = gcc-9 -std=c99 -Wno-unused-result -Wno-unused-function

# how to strip unused code
# STRIPFLAGS = -fdata-sections -ffunction-sections -Wl,--gc-sections -w

# other useful (non-standard) flags
CFLAGS += -g -Wall -pipe

# if you have valgrind, otherwise comment this out
# VALGRIND = valgrind -q --tool=memcheck \
#        --suppressions=$(BASILISK)/openmpi.supp \
#        --leak-check=full

# if gnuplot supports pngcairo, otherwise comment this out
# PNG = pngcairo

# If you have managed to make gdb work (congratulations!), uncomment this
GDB = /usr/local/bin/gdb

# if you don't have md5sum, replace it with something equivalent
GENSUM = shasum
CHECKSUM = shasum -c --status

# OpenGL libraries
# see bview-server.c#installation for explanations
# OPENGLIBS = -lfb_glx -lGLU -lGLEW -lGL -lX11
OPENGLIBS = -L/~/soft/mesa/lib/ -lfb_osmesa -lGLU -lOSMesa

/** You can now compile the libraries provided by basilisk */

cd gl
CFLAGS="-I/Users/antko/soft/mesa/include/ -I/opt/X11/include/ -I/usr/local/include/" make libglutils.a libfb_osmesa.a
cd ..
make -k 
make 
make bview-servers


/** To use bview you'll then need to choose python 2.7 in your terminal : */

conda create --name basilisk python=2.7
conda activate basilisk

/** */

/** */

