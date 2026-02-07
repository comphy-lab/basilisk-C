#!/usr/bin/env python3
"""
Simple VTU inspector: prints coordinates and topology
"""
import struct
import sys
import os

def read_vtu(filename):
    """Read VTU file and return points and topology"""
    with open(filename, 'rb') as f:
        content = f.read()

    # Parse XML metadata
    xml_part = content.split(b'<AppendedData')[0].decode('utf-8')
    pi_idx = xml_part.find('<Piece')
    
    np_idx = xml_part.find('NumberOfPoints="', pi_idx) + 16
    np_end = xml_part.find('"', np_idx)
    num_points = int(xml_part[np_idx:np_end])
    
    nc_idx = xml_part.find('NumberOfCells="', pi_idx) + 15
    nc_end = xml_part.find('"', nc_idx)
    num_cells = int(xml_part[nc_idx:nc_end])

    # Helper to find data offsets
    def get_offset(name, parent_tag=None):
        if parent_tag:
            start = xml_part.find(f'<{parent_tag}>')
            end = xml_part.find(f'</{parent_tag}>')
            section = xml_part[start:end]
        else:
            section = xml_part
        tag_start = section.find(f'Name="{name}"') if name else section.find('<DataArray')
        off_idx = section.find('offset="', tag_start) + 8
        off_end = section.find('"', off_idx)
        return int(section[off_idx:off_end])

    offset_points = get_offset(None, "Points")
    offset_offsets = get_offset("offsets", "Cells")
    offset_types = get_offset("types", "Cells")
    offset_connectivity = get_offset("connectivity", "Cells")
    
    # Find binary data start
    appended_start = content.find(b'<AppendedData encoding="raw">')
    data_start = content.find(b'_', appended_start) + 1
    
    # Read Points
    point_data_start = data_start + offset_points
    coords_start = point_data_start + 8
    points = []
    for i in range(num_points):
        p_idx = coords_start + i * 24
        x, y, z = struct.unpack('<ddd', content[p_idx:p_idx+24])
        points.append((x, y, z))
    
    # Read Offsets
    off_block_start = data_start + offset_offsets
    off_len = struct.unpack('<Q', content[off_block_start:off_block_start+8])[0]
    offsets = list(struct.unpack(f'<{off_len//8}q', content[off_block_start+8:off_block_start+8+off_len]))
    
    # Read Types
    type_block_start = data_start + offset_types
    type_len = struct.unpack('<Q', content[type_block_start:type_block_start+8])[0]
    cell_types = list(struct.unpack(f'<{type_len}b', content[type_block_start+8:type_block_start+8+type_len]))
    
    # Read Connectivity
    conn_block_start = data_start + offset_connectivity
    conn_len = struct.unpack('<Q', content[conn_block_start:conn_block_start+8])[0]
    connectivity = list(struct.unpack(f'<{conn_len//8}q', content[conn_block_start+8:conn_block_start+8+conn_len]))

    return {
        'num_points': num_points,
        'num_cells': num_cells,
        'points': points,
        'offsets': offsets,
        'types': cell_types,
        'connectivity': connectivity
    }

def print_vtu_info(filename):
    """Print coordinates and topology for a VTU file"""
    print(f"\n{'='*70}")
    print(f"File: {os.path.basename(filename)}")
    print(f"{'='*70}")
    
    data = read_vtu(filename)
    
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
        start_ptr = 0 if c_idx == 0 else data['offsets'][c_idx-1]
        end_ptr = data['offsets'][c_idx]
        point_indices = data['connectivity'][start_ptr:end_ptr]
        
        # Check for invalid indices
        invalid = [pid for pid in point_indices if pid < 0 or pid >= data['num_points']]
        if invalid:
            indices_str = f"{point_indices} ⚠ INVALID: {invalid}"
        else:
            indices_str = str(point_indices)
        
        print(f"{c_idx:<6} | {indices_str:<30} | {data['types'][c_idx]:<6}")
    
    print(f"{'='*70}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 inspect_vtu_points.py <file1.vtu> [file2.vtu ...]")
        sys.exit(1)
    
    for filename in sys.argv[1:]:
        if os.path.exists(filename):
            print_vtu_info(filename)
        else:
            print(f"Error: {filename} not found")
