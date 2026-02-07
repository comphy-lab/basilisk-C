#!/usr/bin/env python3
"""
Simple VTKHDF inspector: prints coordinates and topology
"""
import h5py
import sys
import os

def read_vtkhdf(filename):
    """Read VTKHDF file and return points and topology"""
    with h5py.File(filename, 'r') as h5f:
        # VTKHDF format stores data in the VTKHDF group
        if 'VTKHDF' not in h5f:
            raise ValueError(f"File {filename} is not a valid VTKHDF file (missing VTKHDF group)")
        
        ghdf = h5f['VTKHDF']
        
        # Read points (Points)
        points_data = ghdf['Points'][:]
        num_points = len(points_data)
        points = [(p[0], p[1], p[2]) if len(p) == 3 else (p[0], p[1], 0.0) for p in points_data]
        
        # Read connectivity
        connectivity = ghdf['Connectivity'][:].flatten()
        num_connectivity = len(connectivity)
        
        # Read offsets
        offsets = ghdf['Offsets'][:].flatten()
        num_cells = len(offsets) - 1
        
        # Read types
        types = ghdf['Types'][:].flatten()
        
    return {
        'num_points': num_points,
        'num_cells': num_cells,
        'points': points,
        'connectivity': connectivity,
        'offsets': offsets,
        'types': types
    }

def print_vtkhdf_info(filename):
    """Print coordinates and topology for a VTKHDF file"""
    print(f"\n{'='*70}")
    print(f"File: {os.path.basename(filename)}")
    print(f"{'='*70}")
    
    try:
        data = read_vtkhdf(filename)
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    # Print coordinates
    print(f"\nCOORDINATES ({data['num_points']} points):")
    print(f"{'Index':<6} | {'X':<15} | {'Y':<15} | {'Z':<15}")
    print("-" * 55)
    for i, (x, y, z) in enumerate(data['points']):
        print(f"{i:<6} | {x:<15.8f} | {y:<15.8f} | {z:<15.8f}")
    
    # Print topology
    print(f"\nTOPOLOGY ({data['num_cells']} cells):")
    print(f"{'Cell':<6} | {'Point Indices':<30} | {'Type':<6}")
    print("-" * 50)
    for c_idx in range(data['num_cells']):
        start_ptr = data['offsets'][c_idx]
        end_ptr = data['offsets'][c_idx+1]
        point_indices = data['connectivity'][start_ptr:end_ptr]
        
        # Check for invalid indices
        invalid = [pid for pid in point_indices if pid < 0 or pid >= data['num_points']]
        if invalid:
            indices_str = f"{list(point_indices)} ⚠ INVALID: {invalid}"
        else:
            indices_str = str(list(point_indices))
        
        print(f"{c_idx:<6} | {indices_str:<30} | {data['types'][c_idx]:<6}")
    
    print(f"{'='*70}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 test_output_vtkhdf.py <file1.hdf> [file2.hdf ...]")
        sys.exit(1)
    
    for filename in sys.argv[1:]:
        if os.path.exists(filename):
            print_vtkhdf_info(filename)
        else:
            print(f"Error: {filename} not found")
