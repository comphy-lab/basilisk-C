xyz2kdt is a command-line utility used to create the
[Kd-tree-indexed](http://www.wikipedia.org/wiki/Kd-tree) terrain
databases used as input for the objects of the Terrain module. A
summary of the command-line syntax is given by

~~~bash
% xyz2kdt -h
Usage: xyz2kdt [OPTION] BASENAME

Converts the x, y and z coordinates on standard input to a
2D-tree-indexed database suitable for use with the
terrain module of Gerris.

-p N  --pagesize=N  sets the pagesize in bytes (default is 4096)
-v    --verbose     display progress bar
-h    --help        display this help and exit

Report bugs to popinet@users.sf.net
~~~

The format of the data on standard input should look like

    3.501 5.634 -2
    4.601 7.6778 3.456
    ...

where the first field is the value of the *x*-coordinate, the second
field the *y*-coordinate and the last field the *z*-coordinate.
`xyz2kdt` will stop at the first line which does not fit this format.
You may want to check (e.g. using the `--verbose` option) that the
number of points processed matches what you expect. Note that the
database does not enforce any other convention.

# Example 1: building a global terrain topography database using the ETOPO2 dataset

The [ETOPO2 dataset](http://www.ngdc.noaa.gov/mgg/global/)
contains topographic information for the entire surface of the
Earth (both above and below sea level) at a nominal resolution of two
arc-minutes (\~4 km).

The first step is to get the raw data e.g.

~~~bash
% wget http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2/ETOPO2v2-2006/ETOPO2v2c/raw_binary/ETOPO2v2c_i2_LSB.zip
% unzip ETOPO2v2c_i2_LSB.zip
~~~

This is a binary file with a format described in the `ETOPO2v2c_i2_LSB.hdr` file

~~~bash
% cat ETOPO2v2c_i2_LSB.hdr

NCOLS 10800
NROWS 5400
XLLCORNER -180.000000
YLLCORNER -90.000000
CELLSIZE 0.03333333333
NODATA_VALUE -32768
BYTEORDER LSBFIRST
NUMBERTYPE 2_BYTE_INTEGER
MIN_VALUE -10791
MAX_VALUE 8440
~~~

We need to convert this binary file to a text file. We also want to use
the database together with a cartographic projection defined using the
[Map module](object_hierarchy.html#Map_module "Object hierarchy"). By
definition this means that our *x*-, *y*- and *z*-coordinates need to be
the east-positive longitude, north-positive latitude and elevation in
metres. We can easily get these coordinates in a text format suitable
for input into `xyz2kdt` using the following C code:

~~~c
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <arpa/inet.h>

// check that this matches ETOPO2v2c_i2_LSB.hdr

#define NCOLS 10800
#define NROWS 5400
#define XLLCORNER -180.000000
#define YLLCORNER -90.000000
#define CELLSIZE 0.03333333333
#define NODATA_VALUE -32768
#define BYTEORDER LSBFIRST
#define NUMBERTYPE 2_BYTE_INTEGER
#define MIN_VALUE -10791
#define MAX_VALUE 8440

int main (int argc, char * argv[])
{
  double lat, lon;
  int16_t v;
  int i, j;

  for (j = 0; j < NROWS; j++) {
    lat = YLLCORNER + CELLSIZE*j;
    for (i = 0; i < NCOLS; i++) {
      lon = XLLCORNER + CELLSIZE*i;
      assert (fread (&v, sizeof (int16_t), 1, stdin));
      assert (v >= MIN_VALUE && v <= MAX_VALUE);
      printf ("%.8f %.8f %d\n", lon + CELLSIZE/2., - (lat + CELLSIZE/2.), v);
    }
    fprintf (stderr, "\rRow %d/%d              ", j + 1, NROWS);
  }
  fputc ('\n', stderr);
  return 0;
}
~~~

Just copy and paste this code into a file called e.g. `etopo2xyz.c` and
compile using

~~~bash
% cc etopo2xyz.c -o etopo2xyz
~~~

The following command will then read the binary ETOPO2 file, convert it
to the appropriate text format and generate the final terrain database

~~~bash
% ./etopo2xyz < ETOPO2v2c_i2_LSB.bin | xyz2kdt -v etopo2
~~~

If everything went well you should end up (\~ one hour later) with three
large files

~~~bash
% ls etopo2*
etopo2.kdt etopo2.pts etopo2.sum
~~~

which **together** define the terrain database.

Note also that when using several databases simultaneously within
[GfsRefineTerrain](gfsrefineterrain.html "GfsRefineTerrain"), you need
to choose consistent conventions for all the databases (for example a
common [geodetic
system](http://www.wikipedia.org/wiki/Geodetic_system "w:Geodetic_system"){.extiw}
e.g. [WGS84](http://www.wikipedia.org/wiki/WGS84 "w:WGS84"){.extiw}).
The terrain databases do not know anything about projection systems and
it is up to you to enforce your preferred conventions.

# Example 2: building a global terrain dataset using [SRTM data](http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/)

``` {.c}
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <glib.h>

#define NCOLS 1201
#define NROWS 1201
#define CELLSIZE (1./(NCOLS - 1))
#define NODATA_VALUE -32768

int main (int argc, char * argv[])
{
  double lat, lon;
  gint16 v, vs;
  int i, j;

  double xllcorner = atof (argv[1]);
  double yllcorner = atof (argv[2]);

  for (j = 0; j < NROWS; j++) {
    lat = yllcorner + CELLSIZE*(NROWS - 1 - j);
    for (i = 0; i < NCOLS; i++) {
      lon = xllcorner + CELLSIZE*i;
      assert (fread (&v, sizeof (gint16), 1, stdin));
      vs = GINT16_FROM_BE (v);
      if (vs > 0)
    printf ("%.8f %.8f %d\n", lon, lat, vs);
    }
    //    fprintf (stderr, "\rRow %d/%d              ", j + 1, NROWS);
  }
  //  fputc ('\n', stderr);
  return 0;
}
```

~~~bash
% cc `pkg-config glib-2.0 --cflags --libs` -Wall -O2 srtm2kdt.c -o srtm2kdt
~~~

``` {.bash}
for f in dds.cr.usgs.gov/srtm/version2_1/SRTM3/Eurasia/*.hgt.zip; do
    p=`echo $f | sed 's/.*\([NS]\)\([0-9][0-9]\)\([WE]\)\([0-9][0-9][0-9]\)\.hgt\.zip/\1 \2 \3 \4/'`
    lat=`echo $p | awk '{if ($1 == "S") print -$2; else print $2;}'`
    lon=`echo $p | awk '{if ($3 == "W") print -$4; else print $4;}'`
    echo $f >> /dev/stderr
    unzip -c -q $f | ./srtm2kdt $lon $lat
done | xyz2kdt -v srtm_eurasia
```

# Example 3: building a terrain dataset using a GIS raster DEM

You will need [GDAL](http://www.gdal.org). On
debian-like systems this can be installed using

~~~bash
sudo apt-get install gdal-bin
~~~

You can then do something like:

~~~bash
gdal_translate -of XYZ -a_srs EPSG:4326 dem.tif dem.xyz
xyz2kdt -v dem < dem.xyz
rm dem.xyz
~~~

where `EPSG:4326` is
[WGS84](http://www.wikipedia.org/wiki/WGS84). GDAL
supports many formats other than GeoTIFF.
