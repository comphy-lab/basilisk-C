/**
 echo "@include <hdf5.h>" > $BASILISK/ast/std/hdf5.h
 echo "@include <hdf5_hl.h>" > $BASILISK/ast/std/hdf5_hl.h

 and declade the datatypes at the end of $BASILISK/ast/defaults.h

 echo "typedef hid_t, hsize_t, herr_t, H5L_info_t;" >> $BASILISK/ast/defaults.h

 see https://groups.google.com/g/basilisk-fr/c/CM270hBSfWo
 */

#pragma autolink -I${HDF5_INCDIR} -L${HDF5_LIBDIR} -lhdf5

#include <hdf5.h>
#define FAIL -1
#define verbose 1
#define MESG(x) \
	if (verbose)  \
		fprintf(stdout, "%s\n", x);

#define shortcut_slice(n,_alpha)				                \
  double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;\
  if (fabs(alpha) > 0.87)                               \
    continue;                                           \

#define shortcut_periodic()    \
	if (Period.x)                \
		foreach ()                 \
			if (x + Delta > X0 + L0) \
				per_mask[] = 0.;       \
	if (Period.y)                \
		foreach ()                 \
			if (y + Delta > Y0 + L0) \
				per_mask[] = 0.;       \
	if (Period.z)                \
		foreach ()                 \
			if (z + Delta > Z0 + L0) \
				per_mask[] = 0.;


void create_output_h5(char * filename){
	hid_t file_id;
	if (pid() == 0){
		file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		assert(file_id != FAIL);
		MESG("H5Fcreate succeed");
		H5Fclose(file_id);
	}
}

event init(t = 0)
{
	create_output_h5(H5FILE_NAME);
#if dimension == 3
	create_output_h5(H5FILE_NAME2);
#endif
}

trace
void output_xmf(scalar *list, vector *vlist, char * subname, char *stamp)
{
	hid_t acc_tpl1;							/* File access templates */
	hid_t file_id;							/* HDF5 file IDs */
	hid_t group1_id, group2_id; /* HDF5 group IDs */
	hid_t dataspace_id;					/* File dataspace ID */
	hid_t memspace_id;					/* memory dataspace ID */
//hid_t filespace_id;					/* memory dataspace ID */
	hid_t dataset_id;						/* Dataset ID */
	hsize_t dims2[2]; 					/* dataspace dim sizes */
	herr_t status;							/* Generic return value */
//hid_t plist_id;							/* property list identifier */
//H5L_info_t info;
	hsize_t count[2]; /* hyperslab selection parameters */
	hsize_t offset[2];

  /* This a lazy way to match the number of points and vertices when using periodic conditions since vertex are not defined in that case. */
  scalar per_mask[];
  foreach ()
    per_mask[] = 1.;
  shortcut_periodic();

	/* Find the total number of cells/points and per sub-domain */
	int no_point = 0, no_cells = 0, no_point_loc = 0, no_cells_loc = 0;
	foreach_vertex(serial, noauto)
		no_point_loc++;

	foreach (serial, noauto)
		if (per_mask[])
			no_cells_loc++;

	MPI_Allreduce(&no_point_loc, &no_point, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&no_cells_loc, &no_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	/* Get the offsets for each sub-domain */
	int no_point_Loc[npe()], no_point_Off[npe()];
	int no_cells_Loc[npe()], no_cells_Off[npe()];
	for (int i = 0; i < npe(); ++i)
	{
		no_point_Loc[i] = 0.;
		no_cells_Loc[i] = 0.;
	}
	no_point_Loc[pid()] = no_point_loc;
	no_cells_Loc[pid()] = no_cells_loc;

	MPI_Allreduce(no_point_Loc, no_point_Off, npe(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(no_cells_Loc, no_cells_Off, npe(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	offset[0] = 0.;
	if (pid() != 0)
		for (int i = 1; i <= pid(); ++i)
			offset[0] += no_point_Off[i - 1];


	/* Get a marker to reconstruct the topology */
	vertex scalar marker[];
	no_point_loc = 0;
	foreach_vertex(serial, noauto){
#if !TREE
# if dimension == 2
		int _k = (point.i - 2) * ((1 << point.level) + 1) + (point.j - 2);
# else
		int _k = (point.i - 2) * sq((1 << point.level) + 1) + (point.j - 2) * ((1 << point.level) + 1) + (point.k - 2);
# endif
#else // TREE
		int _k = no_point_loc;
#endif
		marker[] = _k + offset[0];
		no_point_loc++;
	}

	/* Write the xmf file. */
	if (pid()==0)
		write_xml_h5_body(list, vlist, H5FILE_NAME, subname, stamp, no_cells, no_point, dimension);


	/* setup file access template with parallel IO access. */
	acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
	assert(acc_tpl1 != FAIL);
	MESG("H5Pcreate access succeed");

	/* set Parallel access with communicator */
	status = H5Pset_fapl_mpio(acc_tpl1, MPI_COMM_WORLD, MPI_INFO_NULL);
	assert(status != FAIL);
	MESG("H5Pset_fapl_mpio succeed");

	/* Open an existing file collectively. */
	file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, acc_tpl1);
	assert(file_id != FAIL);
	MESG("H5Fopen access succeed");

	/* Release file-access template */
	status = H5Pclose(acc_tpl1);
	assert(status != FAIL);

	/* Create (or open) a group for this snapshot. */
	status = H5Lexists(file_id, stamp, H5P_DEFAULT);
	group1_id = (int)status == 0 ? H5Gcreate(file_id, stamp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) : H5Gopen(file_id, stamp, H5P_DEFAULT);
	assert(group1_id != FAIL);
	MESG("H5Gcreate/H5Gopen succeed");

	/****************************************************************************/
	/* Create the dataspace and write the Topology dataset. */
	/* a) Create the dataspace for the dataset. */
	char substamp[1024];
	sprintf(substamp, "/%s/Topology", stamp);
	dims2[0] = no_cells;
	dims2[1] = pow(2, dimension);
	dataspace_id = H5Screate_simple(2, dims2, NULL);
	assert(dataspace_id != FAIL);
	MESG("H5Screate_simple succeed");

	/* b) Create the dataset with default properties and close dataspace. */
	status = H5Lexists(group1_id, substamp, H5P_DEFAULT);
	if ((int)status != 0)
		H5Ldelete(group1_id, substamp, H5P_DEFAULT);
	assert(status != FAIL);
	MESG("H5Lexists succeed");

	dataset_id = H5Dcreate2(group1_id, substamp, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dataset_id != FAIL);
	MESG("H5Dcreate2 succeed");
	H5Sclose(dataspace_id);

	/* c)  Each process defines dataset in memory and writes it to the hyperslab
	in the file. */
	count[0] = no_cells_loc;
	count[1] = pow(2, dimension);
	offset[0] = 0;
	offset[1] = 0;
	if (pid() != 0)
	{
		for (int i = 1; i <= pid(); ++i)
		{
			offset[0] += no_cells_Off[i - 1];
		}
	}
	memspace_id = H5Screate_simple(2, count, NULL);
	assert(memspace_id != FAIL);
	MESG("H5Screate_simple succeed");

	/* d)  Select hyperslab in the file. */
	dataspace_id = H5Dget_space(dataset_id);
	assert(dataspace_id != FAIL);
	MESG("H5Dget_space succeed");

	status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	assert(status != FAIL);
	MESG("H5Sselect_hyperslab succeed");

	long *topo_dset = (long *)malloc(no_cells_loc * count[1] * sizeof(long));
	foreach (serial, noauto){
		if (per_mask[]){
#if !TREE
# if dimension == 2
			int _k = (point.i - 2) * ((1 << point.level)) + (point.j - 2);
# else
			int _k = (point.i - 2) * sq((1 << point.level)) + (point.j - 2) * ((1 << point.level)) + (point.k - 2);
# endif
#endif
			int ii = _k * count[1];
			topo_dset[ii + 0] = (long)marker[];
			topo_dset[ii + 1] = (long)marker[1, 0];
			topo_dset[ii + 2] = (long)marker[1, 1];
			topo_dset[ii + 3] = (long)marker[0, 1];
#if dimension == 3
			topo_dset[ii + 4] = (long)marker[0, 0, 1];
			topo_dset[ii + 5] = (long)marker[1, 0, 1];
			topo_dset[ii + 6] = (long)marker[1, 1, 1];
			topo_dset[ii + 7] = (long)marker[0, 1, 1];
#endif
		}
	}


	/* e) Create property list for collective dataset write. */
	acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
	assert(status != FAIL);
	MESG("H5Pset_dxpl_mpio succeed");

	status = H5Dwrite(dataset_id, H5T_NATIVE_LONG, memspace_id, dataspace_id, acc_tpl1, topo_dset);
	assert(status != FAIL);
	MESG("H5Dwrite succeed");

	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Sclose(memspace_id);
	H5Pclose(acc_tpl1);
	/****************************************************************************/

	/****************************************************************************/
	/* Check if subgroup Geometry exists already */
	status = H5Lexists(group1_id, "Geometry", H5P_DEFAULT);
	if ((int)status == 0)
		group2_id = H5Gcreate(group1_id, "Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	else
	{
		group2_id = H5Gopen(group1_id, "Geometry", H5P_DEFAULT);
	}

	/****************************************************************************/
	/* Create the dataspace and write the Points dataset. */
	/* a) Create the dataspace for the dataset. */
	sprintf(substamp, "/%s/Geometry/Points", stamp);
	dims2[0] = no_point;
	dims2[1] = 3;

	dataspace_id = H5Screate_simple(2, dims2, NULL);
	assert(dataspace_id != FAIL);
	MESG("H5Screate_simple succeed");

	/* b) Create the dataset with default properties and close dataspace. */
	status = H5Lexists(group2_id, substamp, H5P_DEFAULT);
	if ((int)status != 0)
		H5Ldelete(group2_id, substamp, H5P_DEFAULT);
	assert(status != FAIL);
	MESG("H5Lexists succeed");
	dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dataset_id != FAIL);
	MESG("H5Dcreate2 succeed");
	H5Sclose(dataspace_id);

	/* c) Each process defines dataset in memory and writes it to the hyperslab
		 in the file. */
	count[0] = no_point_loc;
	count[1] = 3;
	offset[0] = 0;
	offset[1] = 0;
	if (pid() != 0)
	{
		for (int i = 1; i <= pid(); ++i)
		{
			offset[0] += no_point_Off[i - 1];
		}
	}

	memspace_id = H5Screate_simple(2, count, NULL);
	assert(memspace_id != FAIL);
	MESG("H5Screate_simple succeed");

	/* d)  Select hyperslab in the file. */
	dataspace_id = H5Dget_space(dataset_id);
	assert(dataspace_id != FAIL);
	MESG("H5Dget_space succeed");

	status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	assert(status != FAIL);
	MESG("H5Sselect_hyperslab succeed");

	double *points_dset = (double *)malloc(no_point_loc * 3 * sizeof(double));
	foreach_vertex(serial, noauto){
#if !TREE
# if dimension == 2
		int _k = (point.i - 2) * ((1 << point.level) + 1) + (point.j - 2);
# else
		int _k = (point.i - 2) * sq((1 << point.level) + 1) + (point.j - 2) * ((1 << point.level) + 1) + (point.k - 2);
# endif
#endif
		int ii = _k * 3;
		points_dset[ii + 0] = x;
		points_dset[ii + 1] = y;
#if dimension == 2
		points_dset[ii + 2] = 0.;
#endif
#if dimension == 3
		points_dset[ii + 2] = z;
#endif
	}

	/* e) Create property list for collective dataset write. */
	acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
	assert(status != FAIL);
	MESG("H5Pset_dxpl_mpio succeed");

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, acc_tpl1, points_dset);
	assert(status != FAIL);
	MESG("H5Dwrite succeed");

	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Sclose(memspace_id);
	H5Pclose(acc_tpl1);
	H5Gclose(group2_id);

	/****************************************************************************/
	/* Check if subgroup Cell exists already */
	status = H5Lexists(group1_id, "Cell", H5P_DEFAULT);
	if ((int)status == 0)
		group2_id = H5Gcreate(group1_id, "Cell", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	else
	{
		group2_id = H5Gopen(group1_id, "Cell", H5P_DEFAULT);
	}

	/****************************************************************************/
	/* Create the dataspace and write the scalar dataset. */
	/* Write the scalar datasets. */
	double *scalar_dset = (double *)malloc(no_cells_loc * 1 * sizeof(double));
	for (scalar s in list)
	{
		/* a) Create the dataspace for the dataset. */
		sprintf(substamp, "/%s/Cell/%s", stamp, s.name);
		dims2[0] = no_cells;
		dims2[1] = 1;

		dataspace_id = H5Screate_simple(2, dims2, NULL);
		assert(dataspace_id != FAIL);
		MESG("H5Screate_simple succeed");

		/* b) Create the dataset with default properties and close dataspace. */
		status = H5Lexists(group2_id, substamp, H5P_DEFAULT);
		if ((int)status != 0)
			H5Ldelete(group2_id, substamp, H5P_DEFAULT);
		assert(status != FAIL);
		MESG("H5Lexists succeed");
		dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		assert(dataset_id != FAIL);
		MESG("H5Dcreate2 succeed");
		H5Sclose(dataspace_id);

		/* c)  Each process defines dataset in memory and writes it to the hyperslab
			 in the file. */
		count[0] = no_cells_loc;
		count[1] = 1;
		offset[0] = 0;
		offset[1] = 0;
		if (pid() != 0)
		{
			for (int i = 1; i <= pid(); ++i)
			{
				offset[0] += no_cells_Off[i - 1];
			}
		}
		memspace_id = H5Screate_simple(2, count, NULL);
		assert(memspace_id != FAIL);
		MESG("H5Screate_simple succeed");

		/* d)  Select hyperslab in the file. */
		dataspace_id = H5Dget_space(dataset_id);
		assert(dataspace_id != FAIL);
		MESG("H5Dget_space succeed");

		status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
		assert(status != FAIL);
		MESG("H5Sselect_hyperslab succeed");

		foreach (serial, noauto){
			if (per_mask[]) {
#if !TREE
# if dimension == 2
				int _k = (point.i - 2) * ((1 << point.level)) + (point.j - 2);
# else
				int _k = (point.i - 2) * sq((1 << point.level)) + (point.j - 2) * ((1 << point.level)) + (point.k - 2);
# endif
#endif
				scalar_dset[_k] = s[];
			}
		}

		/* e) Create property list for collective dataset write. */
		acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
		status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
		assert(status != FAIL);
		MESG("H5Pset_dxpl_mpio succeed");

		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, acc_tpl1, scalar_dset);
		assert(status != FAIL);
		MESG("H5Dwrite succeed");

		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		H5Sclose(memspace_id);
		H5Pclose(acc_tpl1);
	}
	/****************************************************************************/


	/****************************************************************************/
	/* Create the dataspace and write the vector datasets. */
	/* Write the scalar datasets. */
	double *vector_dset = (double *)malloc(no_cells_loc * 3 * sizeof(double));
	for (vector v in vlist)
	{
		/* a) Create the dataspace for the dataset. */
		sprintf(substamp, "/%s/Cell/%s", stamp, v.x.name);
		dims2[0] = no_cells;
		dims2[1] = 3;

		dataspace_id = H5Screate_simple(2, dims2, NULL);
		assert(dataspace_id != FAIL);
		MESG("H5Screate_simple succeed");

		/* b) Create the dataset with default properties and close dataspace. */
		status = H5Lexists(group2_id, substamp, H5P_DEFAULT);
		if ((int)status != 0)
			H5Ldelete(group2_id, substamp, H5P_DEFAULT);
		assert(status != FAIL);
		MESG("H5Lexists succeed");
		dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		assert(dataset_id != FAIL);
		MESG("H5Dcreate2 succeed");
		H5Sclose(dataspace_id);

		/* c)  Each process defines dataset in memory and writes it to the hyperslab
			 in the file. */

		count[0] = no_cells_loc;
		count[1] = 3;
		offset[0] = 0;
		offset[1] = 0;
		if (pid() != 0)
		{
			for (int i = 1; i <= pid(); ++i)
			{
				offset[0] += no_cells_Off[i - 1];
			}
		}

		memspace_id = H5Screate_simple(2, count, NULL);
		assert(memspace_id != FAIL);
		MESG("H5Screate_simple succeed");

		/* d)  Select hyperslab in the file. */
		dataspace_id = H5Dget_space(dataset_id);
		assert(dataspace_id != FAIL);
		MESG("H5Dget_space succeed");

		status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
		assert(status != FAIL);
		MESG("H5Sselect_hyperslab succeed");

		foreach (serial, noauto) {
			if (per_mask[]) {
#if !TREE
# if dimension == 2
				int _k = (point.i - 2) * ((1 << point.level)) + (point.j - 2);
# else
				int _k = (point.i - 2) * sq((1 << point.level)) + (point.j - 2) * ((1 << point.level)) + (point.k - 2);
# endif
#endif
				int ii = _k*3;
				vector_dset[ii + 0] = v.x[];
				vector_dset[ii + 1] = v.y[];
#if dimension == 2
				vector_dset[ii + 2] = 0.;
#else
				vector_dset[ii + 2] = v.z[];
#endif
			}
		}

		/* e) Create property list for collective dataset write. */
		acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
		status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
		assert(status != FAIL);
		MESG("H5Pset_dxpl_mpio succeed");

		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, acc_tpl1, vector_dset);
		assert(status != FAIL);
		MESG("H5Dwrite succeed");

		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		H5Sclose(memspace_id);
		H5Pclose(acc_tpl1);
	}
	/****************************************************************************/
	H5Gclose(group2_id);
	H5Gclose(group1_id);
	H5Fclose(file_id);

	free(topo_dset);
	free(points_dset);
	free(scalar_dset);
	free(vector_dset);

	MPI_Barrier(MPI_COMM_WORLD);
}







#if dimension > 2
void output_slice_xmf(scalar *list, vector *vlist, char * subname, char *stamp, coord n, double _alpha)
{
	hid_t acc_tpl1;							/* File access templates */
	hid_t file_id;							/* HDF5 file IDs */
	hid_t group1_id, group2_id; /* HDF5 group IDs */
	hid_t dataspace_id;					/* File dataspace ID */
	hid_t memspace_id;					/* memory dataspace ID */
//hid_t filespace_id;					/* memory dataspace ID */
	hid_t dataset_id;						/* Dataset ID */
	hsize_t dims2[2]; 					/* dataspace dim sizes */
	herr_t status;							/* Generic return value */
//hid_t plist_id;							/* property list identifier */
//H5L_info_t info;
	hsize_t count[2]; /* hyperslab selection parameters */
	hsize_t offset[2];

  /* This a lazy way to match the number of points and vertices when using periodic conditions since vertex are not defined in that case. */
  scalar per_mask[];
	foreach(){
    per_mask[] = 0.;
    shortcut_slice(n,_alpha);
    if (alpha > 0.)
      per_mask[] = 1.;
  }

  shortcut_periodic();

	/* Find the total number of cells/points and per sub-domain */
	vertex scalar marker[];
	int no_point = 0, no_cells = 0, no_point_loc = 0, no_cells_loc = 0;
	foreach_vertex(serial, noauto){
    marker[] = 0.;
    shortcut_slice(n,_alpha);
    marker[] = no_point_loc;
    no_point_loc += 1;
	}

	foreach (serial, noauto)
		if (per_mask[])
			no_cells_loc += 1;

	MPI_Allreduce(&no_point_loc, &no_point, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&no_cells_loc, &no_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	/* Get the offsets for each sub-domain */
	int no_point_Loc[npe()], no_point_Off[npe()];
	int no_cells_Loc[npe()], no_cells_Off[npe()];
	for (int i = 0; i < npe(); ++i)
	{
		no_point_Loc[i] = 0.;
		no_cells_Loc[i] = 0.;
	}
	no_point_Loc[pid()] = no_point_loc;
	no_cells_Loc[pid()] = no_cells_loc;

	MPI_Allreduce(no_point_Loc, no_point_Off, npe(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(no_cells_Loc, no_cells_Off, npe(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	offset[0] = 0.;
	if (pid() != 0)
	{
		for (int i = 1; i <= pid(); ++i)
		{
			offset[0] += no_point_Off[i - 1];
		}
	}

	/* Get a marker to reconstruct the topology */
	foreach_vertex()
		marker[] += offset[0];

	/* Write the xmf file. */
	if (pid()==0)
		write_xml_h5_body(list, vlist, H5FILE_NAME2, subname, stamp, no_cells, no_point, 2);

	/* setup file access template with parallel IO access. */
	acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
	assert(acc_tpl1 != FAIL);
	MESG("H5Pcreate access succeed");

	/* set Parallel access with communicator */
	status = H5Pset_fapl_mpio(acc_tpl1, MPI_COMM_WORLD, MPI_INFO_NULL);
	assert(status != FAIL);
	MESG("H5Pset_fapl_mpio succeed");

	/* Open an existing file collectively. */
	file_id = H5Fopen(H5FILE_NAME2, H5F_ACC_RDWR, acc_tpl1);
	assert(file_id != FAIL);
	MESG("H5Fopen access succeed");

	/* Release file-access template */
	status = H5Pclose(acc_tpl1);
	assert(status != FAIL);

	/* Create (or open) a group for this snapshot. */
	status = H5Lexists(file_id, stamp, H5P_DEFAULT);
	group1_id = (int)status == 0 ? H5Gcreate(file_id, stamp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) : H5Gopen(file_id, stamp, H5P_DEFAULT);
	assert(group1_id != FAIL);
	MESG("H5Gcreate/H5Gopen succeed");

	/****************************************************************************/
	/* Create the dataspace and write the Topology dataset. */
	/* a) Create the dataspace for the dataset. */
	char substamp[1024];
	sprintf(substamp, "/%s/Topology", stamp);
	dims2[0] = no_cells;
	dims2[1] = 4;
	dataspace_id = H5Screate_simple(2, dims2, NULL);
	assert(dataspace_id != FAIL);
	MESG("H5Screate_simple succeed");

	/* b) Create the dataset with default properties and close dataspace. */
	status = H5Lexists(group1_id, substamp, H5P_DEFAULT);
	if ((int)status != 0)
		H5Ldelete(group1_id, substamp, H5P_DEFAULT);
	assert(status != FAIL);
	MESG("H5Lexists succeed");

	dataset_id = H5Dcreate2(group1_id, substamp, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dataset_id != FAIL);
	MESG("H5Dcreate2 succeed");
	H5Sclose(dataspace_id);

	/* c)  Each process defines dataset in memory and writes it to the hyperslab
	in the file. */
	count[0] = no_cells_loc;
	count[1] = 4;
	offset[0] = 0;
	offset[1] = 0;
	if (pid() != 0)
	{
		for (int i = 1; i <= pid(); ++i)
		{
			offset[0] += no_cells_Off[i - 1];
		}
	}
	memspace_id = H5Screate_simple(2, count, NULL);
	assert(memspace_id != FAIL);
	MESG("H5Screate_simple succeed");

	/* d)  Select hyperslab in the file. */
	dataspace_id = H5Dget_space(dataset_id);
	assert(dataspace_id != FAIL);
	MESG("H5Dget_space succeed");

	status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	assert(status != FAIL);
	MESG("H5Sselect_hyperslab succeed");

	long *topo_dset = (long *)malloc(no_cells_loc * count[1] * sizeof(long));
	no_cells_loc = 0;
	foreach (serial, noauto){
		if (per_mask[]) {
			int ii = no_cells_loc * count[1];
			if (n.x == 1) {
				topo_dset[ii + 0] = (long)marker[1, 0, 0];
				topo_dset[ii + 1] = (long)marker[1, 1, 0];
				topo_dset[ii + 2] = (long)marker[1, 1, 1];
				topo_dset[ii + 3] = (long)marker[1, 0, 1];
			}
			else if (n.y == 1) {
				topo_dset[ii + 0] = (long)marker[0, 1, 0];
				topo_dset[ii + 1] = (long)marker[1, 1, 0];
				topo_dset[ii + 2] = (long)marker[1, 1, 1];
				topo_dset[ii + 3] = (long)marker[0, 1, 1];
			}
			else {
				topo_dset[ii + 0] = (long)marker[0, 0, 1];
				topo_dset[ii + 1] = (long)marker[1, 0, 1];
				topo_dset[ii + 2] = (long)marker[1, 1, 1];
				topo_dset[ii + 3] = (long)marker[0, 1, 1];
			}
			no_cells_loc++;
		}
	}

	/* e) Create property list for collective dataset write. */
	acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
	assert(status != FAIL);
	MESG("H5Pset_dxpl_mpio succeed");

	status = H5Dwrite(dataset_id, H5T_NATIVE_LONG, memspace_id, dataspace_id, acc_tpl1, topo_dset);
	assert(status != FAIL);
	MESG("H5Dwrite succeed");

	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Sclose(memspace_id);
	H5Pclose(acc_tpl1);
	/****************************************************************************/

	/****************************************************************************/
	/* Check if subgroup Geometry exists already */
	status = H5Lexists(group1_id, "Geometry", H5P_DEFAULT);
	if ((int)status == 0)
		group2_id = H5Gcreate(group1_id, "Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	else
	{
		group2_id = H5Gopen(group1_id, "Geometry", H5P_DEFAULT);
	}

	/****************************************************************************/
	/* Create the dataspace and write the Points dataset. */
	/* a) Create the dataspace for the dataset. */
	sprintf(substamp, "/%s/Geometry/Points", stamp);
	dims2[0] = no_point;
	dims2[1] = 3;

	dataspace_id = H5Screate_simple(2, dims2, NULL);
	assert(dataspace_id != FAIL);
	MESG("H5Screate_simple succeed");

	/* b) Create the dataset with default properties and close dataspace. */
	status = H5Lexists(group2_id, substamp, H5P_DEFAULT);
	if ((int)status != 0)
		H5Ldelete(group2_id, substamp, H5P_DEFAULT);
	assert(status != FAIL);
	MESG("H5Lexists succeed");
	dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dataset_id != FAIL);
	MESG("H5Dcreate2 succeed");
	H5Sclose(dataspace_id);

	/* c) Each process defines dataset in memory and writes it to the hyperslab
		 in the file. */
	count[0] = no_point_loc;
	count[1] = 3;
	offset[0] = 0;
	offset[1] = 0;
	if (pid() != 0)
	{
		for (int i = 1; i <= pid(); ++i)
		{
			offset[0] += no_point_Off[i - 1];
		}
	}

	memspace_id = H5Screate_simple(2, count, NULL);
	assert(memspace_id != FAIL);
	MESG("H5Screate_simple succeed");

	/* d)  Select hyperslab in the file. */
	dataspace_id = H5Dget_space(dataset_id);
	assert(dataspace_id != FAIL);
	MESG("H5Dget_space succeed");

	status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	assert(status != FAIL);
	MESG("H5Sselect_hyperslab succeed");

	double *points_dset = (double *)malloc(no_point_loc * 3 * sizeof(double));
	no_point_loc = 0;
	foreach_vertex(serial, noauto){
		shortcut_slice(n,_alpha);
		int ii = no_point_loc * 3;
		points_dset[ii + 0] = x;
		points_dset[ii + 1] = y;
		points_dset[ii + 2] = z;
		no_point_loc++;
	}

	/* e) Create property list for collective dataset write. */
	acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
	assert(status != FAIL);
	MESG("H5Pset_dxpl_mpio succeed");

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, acc_tpl1, points_dset);
	assert(status != FAIL);
	MESG("H5Dwrite succeed");

	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Sclose(memspace_id);
	H5Pclose(acc_tpl1);
	H5Gclose(group2_id);

	/****************************************************************************/
	/* Check if subgroup Cell exists already */
	status = H5Lexists(group1_id, "Cell", H5P_DEFAULT);
	if ((int)status == 0)
		group2_id = H5Gcreate(group1_id, "Cell", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	else
	{
		group2_id = H5Gopen(group1_id, "Cell", H5P_DEFAULT);
	}

	/****************************************************************************/
	/* Create the dataspace and write the scalar dataset. */
	/* Write the scalar datasets. */
	double *scalar_dset = (double *)malloc(no_cells_loc * 1 * sizeof(double));
	for (scalar s in list)
	{
		/* a) Create the dataspace for the dataset. */
		sprintf(substamp, "/%s/Cell/%s", stamp, s.name);
		dims2[0] = no_cells;
		dims2[1] = 1;

		dataspace_id = H5Screate_simple(2, dims2, NULL);
		assert(dataspace_id != FAIL);
		MESG("H5Screate_simple succeed");

		/* b) Create the dataset with default properties and close dataspace. */
		status = H5Lexists(group2_id, substamp, H5P_DEFAULT);
		if ((int)status != 0)
			H5Ldelete(group2_id, substamp, H5P_DEFAULT);
		assert(status != FAIL);
		MESG("H5Lexists succeed");
		dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		assert(dataset_id != FAIL);
		MESG("H5Dcreate2 succeed");
		H5Sclose(dataspace_id);

		/* c)  Each process defines dataset in memory and writes it to the hyperslab
			 in the file. */
		count[0] = no_cells_loc;
		count[1] = 1;
		offset[0] = 0;
		offset[1] = 0;
		if (pid() != 0)
		{
			for (int i = 1; i <= pid(); ++i)
			{
				offset[0] += no_cells_Off[i - 1];
			}
		}
		memspace_id = H5Screate_simple(2, count, NULL);
		assert(memspace_id != FAIL);
		MESG("H5Screate_simple succeed");

		/* d)  Select hyperslab in the file. */
		dataspace_id = H5Dget_space(dataset_id);
		assert(dataspace_id != FAIL);
		MESG("H5Dget_space succeed");

		status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
		assert(status != FAIL);
		MESG("H5Sselect_hyperslab succeed");

		no_cells_loc = 0;
		foreach (serial, noauto){
			if (per_mask[]) {
				if (n.x == 1)
					scalar_dset[no_cells_loc] = 0.5*(val(s) + val(s,1,0,0));
				else if (n.y == 1)
					scalar_dset[no_cells_loc] = 0.5*(val(s) + val(s,0,1,0));
				else
					scalar_dset[no_cells_loc] = 0.5*(val(s) + val(s,0,0,1));
				no_cells_loc++;
			}
		}

		/* e) Create property list for collective dataset write. */
		acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
		status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
		assert(status != FAIL);
		MESG("H5Pset_dxpl_mpio succeed");

		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, acc_tpl1, scalar_dset);
		assert(status != FAIL);
		MESG("H5Dwrite succeed");

		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		H5Sclose(memspace_id);
		H5Pclose(acc_tpl1);
	}
	/****************************************************************************/


	/****************************************************************************/
	/* Create the dataspace and write the vector datasets. */
	/* Write the scalar datasets. */
	double *vector_dset = (double *)malloc(no_cells_loc * 3 * sizeof(double));
	for (vector v in vlist)
	{
		/* a) Create the dataspace for the dataset. */
		sprintf(substamp, "/%s/Cell/%s", stamp, v.x.name);
		dims2[0] = no_cells;
		dims2[1] = 3;

		dataspace_id = H5Screate_simple(2, dims2, NULL);
		assert(dataspace_id != FAIL);
		MESG("H5Screate_simple succeed");

		/* b) Create the dataset with default properties and close dataspace. */
		status = H5Lexists(group2_id, substamp, H5P_DEFAULT);
		if ((int)status != 0)
			H5Ldelete(group2_id, substamp, H5P_DEFAULT);
		assert(status != FAIL);
		MESG("H5Lexists succeed");
		dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		assert(dataset_id != FAIL);
		MESG("H5Dcreate2 succeed");
		H5Sclose(dataspace_id);

		/* c)  Each process defines dataset in memory and writes it to the hyperslab
			 in the file. */

		count[0] = no_cells_loc;
		count[1] = 3;
		offset[0] = 0;
		offset[1] = 0;
		if (pid() != 0)
		{
			for (int i = 1; i <= pid(); ++i)
			{
				offset[0] += no_cells_Off[i - 1];
			}
		}

		memspace_id = H5Screate_simple(2, count, NULL);
		assert(memspace_id != FAIL);
		MESG("H5Screate_simple succeed");

		/* d)  Select hyperslab in the file. */
		dataspace_id = H5Dget_space(dataset_id);
		assert(dataspace_id != FAIL);
		MESG("H5Dget_space succeed");

		status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
		assert(status != FAIL);
		MESG("H5Sselect_hyperslab succeed");

		no_cells_loc = 0;
		foreach (serial, noauto){
			if (per_mask[]) {
				int ii = no_cells_loc*3;
				if (n.x == 1) {
					vector_dset[ii + 0] = 0.5*(val(v.x) + val(v.x,1,0,0));
					vector_dset[ii + 1] = 0.5*(val(v.y) + val(v.y,1,0,0));
					vector_dset[ii + 2] = 0.5*(val(v.z) + val(v.z,1,0,0));
				}
				else if (n.y == 1) {
					vector_dset[ii + 0] = 0.5*(val(v.x) + val(v.x,0,1,0));
					vector_dset[ii + 1] = 0.5*(val(v.y) + val(v.y,0,1,0));
					vector_dset[ii + 2] = 0.5*(val(v.z) + val(v.z,0,1,0));
				}
				else {
					vector_dset[ii + 0] = 0.5*(val(v.x) + val(v.x,0,0,1));
					vector_dset[ii + 1] = 0.5*(val(v.y) + val(v.y,0,0,1));
					vector_dset[ii + 2] = 0.5*(val(v.z) + val(v.z,0,0,1));
				}
				no_cells_loc++;
			}
		}

		/* e) Create property list for collective dataset write. */
		acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
		status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
		assert(status != FAIL);
		MESG("H5Pset_dxpl_mpio succeed");

		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, acc_tpl1, vector_dset);
		assert(status != FAIL);
		MESG("H5Dwrite succeed");

		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		H5Sclose(memspace_id);
		H5Pclose(acc_tpl1);
	}
	/****************************************************************************/
	H5Gclose(group2_id);
	H5Gclose(group1_id);
	H5Fclose(file_id);

	free(topo_dset);
	free(points_dset);
	free(scalar_dset);
	free(vector_dset);

	MPI_Barrier(MPI_COMM_WORLD);
}
#endif

#undef shortcut_periodic
#undef shortcut_slice