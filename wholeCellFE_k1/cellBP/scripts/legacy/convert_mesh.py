# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 10:47:08 2023

@author: andre 
"""

with open('disk_20um_diam_gradient.ply2', 'r') as f_in, open('disk_20um_diam_gradient.txt', 'w') as f_out:
    # read and write header lines up to node data
    header1 = f_in.readline().strip()
    
    # write modified header for node data to output file
    f_out.write('Nodes\n')
    f_out.write(header1 + '\n')

    # read number of nodes
    num_nodes = int(header1.split()[-1])

    # read header line up to element data
    header2 = f_in.readline().strip()

    

    # process node data
    for i in range(num_nodes):
        line = f_in.readline()
        # split line into columns, extract x and y coordinates
        columns = line.split()
        x, y, z = columns[:3]
        f_out.write(f'{x} {y}\n')  # write modified line to output file

    # write header for element data to output file
    f_out.write('Elements\n')
    f_out.write(header2 + '\n')

    # read number of elements
    num_elems = int(header2.split()[-1])

    # process element data
    for i in range(num_elems):
        line = f_in.readline()
        # split line into columns, extract node numbers
        columns = line.split()
        elem_type, n1, n2, n3 = columns[:4]
        f_out.write(f'{n1} {n2} {n3}\n')  # write modified line to output file

