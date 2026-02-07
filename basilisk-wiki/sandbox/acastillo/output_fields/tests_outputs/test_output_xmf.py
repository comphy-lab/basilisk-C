#!/usr/bin/env python3
"""
Simple XDMF/HDF5 inspector: prints coordinates and topology
"""
import h5py
import xml.etree.ElementTree as ET
import sys
import os

def read_xdmf(xmf_filename):
    """Read XDMF file and return points and topology"""
    
    # Parse the XDMF XML file
    tree = ET.parse(xmf_filename)
    root = tree.getroot()
    
    # Find the Grid element
    domain = root.find('.//Domain')
    grid = domain.find('.//Grid')
    
    # Get the HDF5 filename from the XML
    h5_filename = None
    for data_item in grid.findall('.//DataItem'):
        text = data_item.text.strip()
        if '.h5:' in text:
            h5_filename = text.split(':')[0]
            break
    
    if h5_filename is None:
        raise ValueError("Could not find HDF5 filename in XDMF file")
    
    # Make path relative to XDMF file location
    xmf_dir = os.path.dirname(xmf_filename)
    h5_path = os.path.join(xmf_dir, h5_filename) if xmf_dir else h5_filename
    
    # Read data from HDF5 file
    with h5py.File(h5_path, 'r') as h5f:
        # Read points (geometry)
        points_data = h5f['/Geometry/Points'][:]
        num_points = len(points_data)
        points = [(p[0], p[1], p[2]) for p in points_data]
        
        # Read topology (connectivity)
        # This is stored as a 2D array: num_cells x points_per_cell
        topology_data = h5f['/Topology'][:]
        num_cells = len(topology_data)
        points_per_cell = topology_data.shape[1]
        
        # Convert to list of lists for easier handling (convert to plain Python ints)
        connectivity = [[int(x) for x in row] for row in topology_data]
    
    return {
        'num_points': num_points,
        'num_cells': num_cells,
        'points': points,
        'connectivity': connectivity,
        'points_per_cell': points_per_cell,
        'h5_file': h5_filename
    }

def print_xdmf_info(xmf_filename):
    """Print coordinates and topology for an XDMF file"""
    print(f"\n{'='*70}")
    print(f"File: {os.path.basename(xmf_filename)}")
    print(f"{'='*70}")
    
    data = read_xdmf(xmf_filename)
    
    print(f"HDF5 file: {data['h5_file']}")
    print(f"Points per cell: {data['points_per_cell']}")
    
    # Print coordinates
    print(f"\nCOORDINATES ({data['num_points']} points):")
    print(f"{'Index':<6} | {'X':<15} | {'Y':<15} | {'Z':<15}")
    print("-" * 55)
    for i, (x, y, z) in enumerate(data['points']):
        print(f"{i:<6} | {x:<15.8f} | {y:<15.8f} | {z:<15.8f}")
    
    # Print topology
    print(f"\nTOPOLOGY ({data['num_cells']} cells):")
    print(f"{'Cell':<6} | {'Point Indices':<40}")
    print("-" * 50)
    for c_idx in range(data['num_cells']):
        point_indices = data['connectivity'][c_idx]
        
        # Check for invalid indices
        invalid = [pid for pid in point_indices if pid < 0 or pid >= data['num_points']]
        if invalid:
            indices_str = f"{point_indices} ⚠ INVALID: {invalid}"
        else:
            indices_str = str(point_indices)
        
        print(f"{c_idx:<6} | {indices_str:<40}")
    
    print(f"{'='*70}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 test_output_xmf.py <file1.xmf> [file2.xmf ...]")
        sys.exit(1)
    
    for filename in sys.argv[1:]:
        if os.path.exists(filename):
            print_xdmf_info(filename)
        else:
            print(f"Error: {filename} not found")
