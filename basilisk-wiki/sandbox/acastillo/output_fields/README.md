# Output Fields

This includes routines to write fields into different formats:

- `vtu`: These routines are compatible with the VTK XML file format
  for unstructured grids File formats for VTK version 4.2. Heavy data is appended
  at the end of each file in raw binary. If used on MPI, functions write one
  file per process linked together using a `.pvtu` file. 

- `xdmf`: These routines are compatible with the XDMF Model and Format. Heavy
  data is stored inside a `.h5` file using the Hierarchical Data Format HDF5.
  This is a parallel implementation (only) suitable for HPC. Requires linking to
  `-lhdf5`. To use HDF5 we need to "declare" the library in the `qcc` system

    ```
    echo "@include <hdf5.h>" > $BASILISK/ast/std/hdf5.h
    echo "@include <hdf5_hl.h>" > $BASILISK/ast/std/hdf5_hl.h
    ```
    and declare the datatypes at the end of `$BASILISK/ast/defaults.h`

    ```
    echo "typedef hid_t, hsize_t, herr_t, H5L_info_t;" >> $BASILISK/ast/defaults.h
    ```

    see https://groups.google.com/g/basilisk-fr/c/CM270hBSfWo

- `vtkhdf`: These routines are compatible with the VTKHDF File Format for
  unstructured grids. The `VTKHDF` file format is a file format relying on HDF5.
  It is meant to provide good I/O performance as well as robust and flexible
  parallel I/O capabilities. 

    ```
    /                        Group
    /VTKHDF                  Group
    /VTKHDF/CellData         Group
    /VTKHDF/CellData/f       Dataset {1676032, 1}
    /VTKHDF/CellData/p       Dataset {1676032, 1}
    /VTKHDF/CellData/u.x     Dataset {1676032, 3}
    /VTKHDF/Connectivity     Dataset {13408256, 1}
    /VTKHDF/FieldData        Group
    /VTKHDF/NumberOfCells    Dataset {4}
    /VTKHDF/NumberOfConnectivityIds Dataset {4}
    /VTKHDF/NumberOfPoints   Dataset {4}
    /VTKHDF/Offsets          Dataset {1676036, 1}
    /VTKHDF/PointData        Group
    /VTKHDF/Points           Dataset {1772232, 3}
    /VTKHDF/Types            Dataset {1676032, 1}
    ```

The other routines might be outdated.