/**
# XMDF

Light data is stored using eXtensible Markup Language (XML) to describe the data
model and the data format. We're going to output two main types of files. Heavy
files are set to `H5FILE_NAME` for 3D fields and `H5FILE_NAME2` for extracting
slices.
*/


#ifndef H5FILE_NAME
#define H5FILE_NAME "output_fields.h5"
#endif

#ifndef H5FILE_NAME2
#define H5FILE_NAME2 "output_slices.h5"
#endif

/**
## Light Data
 Light files are created by the auxiliar functions `write_xml_h5_head`,
`write_xml_h5_body` and `write_xml_h5_tail`.
*/

void write_xml_h5_head(char * subname){
	if (pid() == 0){
		char name[111];
		sprintf(name, "%s.xmf", subname);
		FILE *fp = fopen(name, "w");

#if defined(_OPENMP)
		int num_omp = omp_get_max_threads();
		omp_set_num_threads(1);
#endif

		/* Write the xmf file. */
		fputs("<?xml version=\"1.0\"?>\n", fp);
		fputs("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n", fp);
		fputs("<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n",
					fp);
		fputs("\t <Domain>\n", fp);
		fputs("\t\t <Grid GridType=\"Collection\" CollectionType=\"Temporal\" >\n", fp);
		fputs("\n", fp);

#if defined(_OPENMP)
		omp_set_num_threads(num_omp);
#endif
		fflush(fp);
		fclose(fp);
	}
}

void write_xml_h5_body(scalar *list, vector *vlist, char * filename, char * subname, char *stamp, int no_cells, int no_point, int dim){
	char name[111];
	sprintf(name, "%s.xmf", subname);
  FILE * fp = fopen(name, "a");

#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

	/* Write the xmf file. */
	fputs("\t\t\t <Grid GridType=\"Uniform\">\n", fp);
	fprintf(fp, "\t\t\t\t <Time Type=\"Single\" Value=\"%g\" />\n", t);

	if (dim == 2)
	{
		fprintf(fp, "\t\t\t\t <Topology TopologyType=\"Quadrilateral\" Dimensions=\"%d\">\n", no_cells);
		fprintf(fp, "\t\t\t\t\t <DataItem Dimensions=\"%d 4\" NumberType=\"Int\" Precision=\"8\" Format=\"HDF\">\n", no_cells);
	}

	if (dim == 3)
	{
		fprintf(fp, "\t\t\t\t <Topology TopologyType=\"Hexahedron\" Dimensions=\"%d\">\n", no_cells);
		fprintf(fp, "\t\t\t\t\t <DataItem Dimensions=\"%d 8\" NumberType=\"Int\" Precision=\"8\" Format=\"HDF\">\n", no_cells);
	}

	fprintf(fp, "\t\t\t\t\t\t %s:/%s/Topology\n", filename, stamp);
	fputs("\t\t\t\t\t </DataItem>\n", fp);
	fputs("\t\t\t\t </Topology>\n", fp);
	fputs("\t\t\t\t <Geometry GeometryType=\"XYZ\">\n", fp);
	fprintf(fp, "\t\t\t\t\t <DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", no_point);
	fprintf(fp, "\t\t\t\t\t\t %s:/%s/Geometry/Points\n", filename, stamp);
	fputs("\t\t\t\t\t </DataItem>\n", fp);
	fputs("\t\t\t\t </Geometry>\n", fp);
	for (scalar s in list)
	{
		fprintf(fp, "\t\t\t\t <Attribute Name=\"%s\" Active=\"1\" AttributeType=\"Scalar\" Center=\"Cell\">\n", s.name);
		fprintf(fp, "\t\t\t\t\t <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", no_cells);
		fprintf(fp, "\t\t\t\t\t\t %s:/%s/Cell/%s\n", filename, stamp, s.name);
		fputs("\t\t\t\t\t </DataItem>\n", fp);
		fputs("\t\t\t\t </Attribute>\n", fp);
	}
	for (vector v in vlist)
	{
		fprintf(fp, "\t\t\t\t <Attribute Name=\"%s\" Active=\"1\" AttributeType=\"Vector\" Center=\"Cell\">\n", v.x.name);
		fprintf(fp, "\t\t\t\t\t <DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", no_cells);
		fprintf(fp, "\t\t\t\t\t\t %s:/%s/Cell/%s\n", filename, stamp, v.x.name);
		fputs("\t\t\t\t\t </DataItem>\n", fp);
		fputs("\t\t\t\t </Attribute>\n", fp);
	}
	fputs("\t\t\t </Grid>\n", fp);
	fputs("\n", fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
	fflush(fp);
	fclose (fp);
}

void write_xml_h5_tail(char *subname)
{
	if (pid() == 0){
		char name[111];
		sprintf(name, "%s.xmf", subname);
		FILE *fp = fopen(name, "a");

#if defined(_OPENMP)
		int num_omp = omp_get_max_threads();
		omp_set_num_threads(1);
#endif

		/* Write the xmf file. */
		fputs("\t\t </Grid>\n", fp);
		fputs("\t </Domain>\n", fp);
		fputs("</Xdmf>\n", fp);

#if defined(_OPENMP)
		omp_set_num_threads(num_omp);
#endif
		fflush(fp);
		fclose(fp);
	}
}

/**
### Typical example

A typical example would look like this

```
  create_output_h5(H5FILE_NAME);
  write_xml_h5_head("test");
  output_xmf (slist, vlist, "test", "3D000");
  write_xml_h5_tail("test");
```
This creates a file "output_fields.h5" which contains
```
>> h5ls -r output_fields.h5
/                        Group
/3D000                   Group
/3D000/Cell              Group
/3D000/Cell/s            Dataset {3050496, 1}
/3D000/Cell/u.x          Dataset {3050496, 3}
...

```
and a file `test.xmf` explaining how to read the data
```
<?xml version="1.0"?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">
	 <Domain>
		 <Grid GridType="Collection" CollectionType="Temporal" >

			 ...

		 </Grid>
	 </Domain>
</Xdmf>

```

*/

/**
## Heavy Data
Heavy data is composed of large datasets stored
using the Hierarchical Data Format HDF5. Additional features include support
for a large number of objects, file compression, a parallel I/O implementation
through the MPI-IO or MPI POSIX drivers.

To use HDF5 we need to "declare" the library in the qcc system
```
 echo "@include <hdf5.h>" > $BASILISK/ast/std/hdf5.h
 echo "@include <hdf5_hl.h>" > $BASILISK/ast/std/hdf5_hl.h
```

and declade the datatypes at the end of $BASILISK/ast/defaults.h
```
 echo "typedef hid_t, hsize_t, herr_t, H5L_info_t;" >> $BASILISK/ast/defaults.h
```
see more about declaractions with qcc in the
[forum](https://groups.google.com/g/basilisk-fr/c/CM270hBSfWo).
*/
#include "output_hdf5.h"
