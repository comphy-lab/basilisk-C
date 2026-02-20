#!/usr/bin/env python3
"""
xdmf_chunk.py
=============
Extract a contiguous slice of cells from a (potentially enormous) domain.h5 / domain.xmf
pair and write a self-contained chunk HDF5 + XMF file that ParaView can open directly.

The domain stores:
  /Topology          (N_cells, 8)  int64   -- hex-8 connectivity (global point ids)
  /Geometry/Points   (N_pts,   3)  float64 -- XYZ coordinates
  /Cells/<field>     (N_cells,  …) float64 -- one dataset per field

Strategy
--------
1. Read only the topology rows [cell_start, cell_end).
2. Collect the unique global point ids referenced by those cells.
3. Build a local-index remapping (global → local).
4. Read only the required rows of /Geometry/Points.
5. Read the field data for the selected cells.
6. Write everything to chunk_<tag>.h5 and chunk_<tag>.xmf.

Usage
-----
    # absolute cell indices
    python xdmf_chunk.py \\
        --input  domain.h5 \\
        --output chunk \\
        --cell-start 0 \\
        --cell-end   1000000

    # percentage of the total cell count (0–100)
    python xdmf_chunk.py \\
        --input  domain.h5 \\
        --output chunk \\
        --pct-start 25 \\
        --pct-end   50

Optional flags
--------------
    --fields    f p u.x        (default: all datasets found under /Cells)
    --time      0.0           (value written into <Time> element of the XMF)
    --pct-start 25.0          (start as % of total cells; overrides --cell-start)
    --pct-end   50.0          (end   as % of total cells; overrides --cell-end)

Requirements
------------
    pip install h5py numpy
"""

import argparse
import os
import sys

import h5py
import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def detect_fields(h5_file: h5py.File) -> list[str]:
    """Return all dataset names found under /Cells."""
    if "Cells" not in h5_file:
        return []
    return list(h5_file["Cells"].keys())


def xmf_scalar(name: str, n_cells: int, hdf5_filename: str) -> str:
    return (
        f'\t\t\t<Attribute Name="{name}" AttributeType="Scalar" Center="Cell">\n'
        f'\t\t\t\t<DataItem Dimensions="{n_cells}" NumberType="Float" Precision="8" Format="HDF">\n'
        f'\t\t\t\t\t{hdf5_filename}:/Cells/{name}\n'
        f'\t\t\t\t</DataItem>\n'
        f'\t\t\t</Attribute>\n'
    )


def xmf_vector(name: str, n_cells: int, n_comp: int, hdf5_filename: str) -> str:
    return (
        f'\t\t\t<Attribute Name="{name}" AttributeType="Vector" Center="Cell">\n'
        f'\t\t\t\t<DataItem Dimensions="{n_cells} {n_comp}" NumberType="Float" Precision="8" Format="HDF">\n'
        f'\t\t\t\t\t{hdf5_filename}:/Cells/{name}\n'
        f'\t\t\t\t</DataItem>\n'
        f'\t\t\t</Attribute>\n'
    )


def write_xmf(
    xmf_path: str,
    hdf5_filename: str,
    n_cells: int,
    n_points: int,
    node_per_cell: int,
    fields_meta: list[tuple],   # list of (name, shape)  shape = (n_cells,) or (n_cells, k)
    time_value: float = 0.0,
    topology_type: str = "Hexahedron",
):
    """Write the lightweight XMF descriptor file."""

    attr_lines = []
    for name, shape in fields_meta:
        if len(shape) == 1 or (len(shape) == 2 and shape[1] == 1):
            attr_lines.append(xmf_scalar(name, n_cells, hdf5_filename))
        else:
            n_comp = shape[1]
            attr_lines.append(xmf_vector(name, n_cells, n_comp, hdf5_filename))

    attrs = "".join(attr_lines)

    xmf = f"""\
<?xml version="1.0"?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [
<!ENTITY HeavyData "{hdf5_filename}:">
]>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">
\t<Domain>
\t\t<Grid Name="Unstructured Grid" GridType="Uniform">
\t\t\t<Time Type="Single" Value="{time_value}" />
\t\t\t<Topology TopologyType="{topology_type}" NumberOfElements="{n_cells}">
\t\t\t\t<DataItem Format="HDF" Dimensions="{n_cells} {node_per_cell}" DataType="Int" Precision="8" >
\t\t\t\t\t&HeavyData;/Topology
\t\t\t\t</DataItem>
\t\t\t</Topology>
\t\t\t<Geometry GeometryType="XYZ">
\t\t\t\t<DataItem Format="HDF" NumberType="Float" Dimensions="{n_points} 3" Precision="8" >
\t\t\t\t\t&HeavyData;/Geometry/Points
\t\t\t\t</DataItem>
\t\t\t</Geometry>
{attrs}\t\t</Grid>
\t</Domain>
</Xdmf>
"""
    with open(xmf_path, "w") as fh:
        fh.write(xmf)
    print(f"  XMF written  : {xmf_path}")


# ---------------------------------------------------------------------------
# Core extraction
# ---------------------------------------------------------------------------

def extract_chunk(
    input_h5: str,
    output_prefix: str,
    cell_start: int,
    cell_end: int,
    fields: list[str] | None = None,
    time_value: float = 0.0,
    verbose: bool = True,
):
    """
    Parameters
    ----------
    input_h5      : path to the source domain.h5
    output_prefix : prefix for output files, e.g. "chunk" → chunk.h5 + chunk.xmf
    cell_start    : first cell index to extract (inclusive)
    cell_end      : last  cell index to extract (exclusive)
    fields        : list of field names under /Cells to copy; None → all fields
    time_value    : simulation time to embed in XMF
    verbose       : print progress messages
    """

    def log(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    output_dir = os.path.dirname(output_prefix) or "."
    tag = os.path.basename(output_prefix)
    h5_out_path  = os.path.join(output_dir, f"{tag}.h5")
    xmf_out_path = os.path.join(output_dir, f"{tag}.xmf")

    with h5py.File(input_h5, "r") as src:

        # ---- discover fields -------------------------------------------------
        all_fields = detect_fields(src)
        if fields is None:
            fields = all_fields
        else:
            # validate
            missing = [f for f in fields if f not in all_fields]
            if missing:
                sys.exit(f"ERROR: fields not found in /Cells: {missing}\n"
                         f"Available: {all_fields}")

        # ---- topology shape --------------------------------------------------
        topo_ds   = src["/Topology"]
        n_cells_total, node_per_cell = topo_ds.shape
        cell_end  = min(cell_end, n_cells_total)   # clamp to file size
        n_cells   = cell_end - cell_start

        if n_cells <= 0:
            sys.exit(f"ERROR: cell range [{cell_start}, {cell_end}) is empty or invalid "
                     f"(dataset has {n_cells_total} cells).")

        log(f"Source        : {input_h5}")
        log(f"Total cells   : {n_cells_total:,}")
        log(f"Chunk range   : [{cell_start:,}, {cell_end:,})  → {n_cells:,} cells")
        log(f"Node/cell     : {node_per_cell}")

        # ---- 1. read connectivity chunk -------------------------------------
        log("Reading topology chunk …")
        topo_chunk = topo_ds[cell_start:cell_end, :]   # shape (n_cells, node_per_cell)

        # ---- 2. unique point indices ----------------------------------------
        log("Finding unique point indices …")
        global_pt_ids = np.unique(topo_chunk)          # sorted 1-D array
        n_points = len(global_pt_ids)
        log(f"Unique points : {n_points:,}")

        # ---- 3. remap connectivity to local indices -------------------------
        log("Remapping connectivity …")
        # Build a lookup: global_id → local_id
        # np.searchsorted works because global_pt_ids is sorted
        local_topo = np.searchsorted(global_pt_ids, topo_chunk).astype(np.int64)

        # ---- 4. read only required point coordinates -----------------------
        log("Reading point coordinates …")
        pts_ds = src["/Geometry/Points"]

        # HDF5 fancy indexing on large datasets can be very slow; use a loop
        # over sorted consecutive runs when possible.  For simplicity we use
        # a single fancy-index call — h5py handles it efficiently for sorted
        # indices.
        points = pts_ds[global_pt_ids, :]   # shape (n_points, 3)

        # ---- 5. read fields -------------------------------------------------
        fields_data = {}
        for fname in fields:
            log(f"Reading field  : /Cells/{fname} …")
            fields_data[fname] = src[f"/Cells/{fname}"][cell_start:cell_end]

        # ---- 6. write output HDF5 ------------------------------------------
        log(f"Writing HDF5   : {h5_out_path}")
        with h5py.File(h5_out_path, "w") as dst:

            # Topology (local indices)
            dst.create_dataset(
                "/Topology",
                data=local_topo,
                dtype=np.int64,
                chunks=(min(1_048_576, n_cells), node_per_cell),
                compression="gzip",
                compression_opts=6,
                shuffle=True,
            )

            # Geometry
            dst.create_dataset(
                "/Geometry/Points",
                data=points,
                dtype=np.float64,
                chunks=(min(1_048_576, n_points), 3),
                compression="gzip",
                compression_opts=6,
                shuffle=True,
            )

            # Fields
            for fname, data in fields_data.items():
                chunks = (min(1_048_576, n_cells),) + data.shape[1:]
                dst.create_dataset(
                    f"/Cells/{fname}",
                    data=data,
                    dtype=np.float64,
                    chunks=chunks,
                    compression="gzip",
                    compression_opts=6,
                    shuffle=True,
                )

        # ---- 7. write XMF --------------------------------------------------
        topology_type_map = {4: "Quadrilateral", 8: "Hexahedron"}
        topology_type = topology_type_map.get(node_per_cell, f"Unknown_{node_per_cell}")

        fields_meta = [(fname, fields_data[fname].shape) for fname in fields]
        hdf5_basename = os.path.basename(h5_out_path)

        write_xmf(
            xmf_path=xmf_out_path,
            hdf5_filename=hdf5_basename,
            n_cells=n_cells,
            n_points=n_points,
            node_per_cell=node_per_cell,
            fields_meta=fields_meta,
            time_value=time_value,
            topology_type=topology_type,
        )

        log("Done.")
        return h5_out_path, xmf_out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Extract a chunk of cells from a domain.h5/domain.xmf pair.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",  "-i", required=True,
                   help="Path to the source HDF5 file (e.g. domain.h5)")
    p.add_argument("--output", "-o", default="chunk",
                   help="Output path prefix (produces <prefix>.h5 and <prefix>.xmf)")
    p.add_argument("--cell-start", type=int, default=0,
                   help="First cell index to extract (inclusive); ignored if --pct-start is set")
    p.add_argument("--cell-end",   type=int, default=1_000_000,
                   help="Last  cell index to extract (exclusive); ignored if --pct-end is set")
    p.add_argument("--pct-start",  type=float, default=None,
                   help="Start as %% of total cells [0, 100]; overrides --cell-start")
    p.add_argument("--pct-end",    type=float, default=None,
                   help="End   as %% of total cells [0, 100]; overrides --cell-end")
    p.add_argument("--fields", nargs="*", default=None,
                   help="Field names under /Cells to copy (default: all)")
    p.add_argument("--time",   type=float, default=0.0,
                   help="Simulation time to embed in the XMF <Time> element")
    p.add_argument("--quiet",  action="store_true",
                   help="Suppress progress messages")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    cell_start = args.cell_start
    cell_end   = args.cell_end

    # Resolve percentages — requires a quick peek at the dataset shape
    if args.pct_start is not None or args.pct_end is not None:
        with h5py.File(args.input, "r") as _f:
            n_cells_total = _f["/Topology"].shape[0]
        if args.pct_start is not None:
            cell_start = int(round(args.pct_start / 100.0 * n_cells_total))
        if args.pct_end is not None:
            cell_end   = int(round(args.pct_end   / 100.0 * n_cells_total))
        if not args.quiet:
            print(f"Total cells   : {n_cells_total:,}")
            print(f"Resolved range: [{cell_start:,}, {cell_end:,})")

    extract_chunk(
        input_h5=args.input,
        output_prefix=args.output,
        cell_start=cell_start,
        cell_end=cell_end,
        fields=args.fields if args.fields else None,
        time_value=args.time,
        verbose=not args.quiet,
    )
