#!/usr/bin/python
import numpy as np
from tvtk.api import tvtk
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter


if __name__ == '__main__':
    parser = ArgumentParser(description='fix vtk mesh in doubles for FreeSurfer',
                            formatter_class=RawTextHelpFormatter)
    g_input = parser.add_argument_group('Inputs')
    g_input.add_argument('-i','--input', action='store', help='input file')

    g_output = parser.add_argument_group('Outputs')
    g_output.add_argument('-o', '--output', action='store', help='output file')
    opts = parser.parse_args()

    r = tvtk.PolyDataReader(file_name=opts.input)
    r.update()
    m = r.output
    
    points = np.array(m.points).astype(np.float32)
    m.points = points
    
    w = tvtk.PolyDataWriter(file_name=opts.output)
    w.set_input_data(m)
    w.write()
