import numpy as np
import nibabel as nib
from tvtk.api import tvtk
import datetime as dt
import os
import argparse
import sys
import subprocess as sp


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!"%arg)
    else:
        return arg

def transform(origin, target, tfm, output_prefix, use_paraview=False ):
    del_vtk = False
    reg = np.genfromtxt(tfm,skiprows=4,comments='r')
    cmd_info = "mri_info --ras2vox-tkr %s" % target
    proc = sp.Popen( cmd_info, stdout=sp.PIPE, shell=True )
    data = bytearray( proc.stdout.read() )
    dti_ras2vox_tkr = np.reshape( np.fromstring(data.decode("utf-8"), sep='\n' ), (4,-1) )
    ac = nib.load( target ).get_affine() # this is the same as dti_vox2ras
#    tkr_to_b0 = np.dot( np.linalg.inv(dti_ras2vox_tkr), ac )
    tkr_to_b0 = np.dot( ac, dti_ras2vox_tkr )
#    xfm = np.dot( reg, tkr_to_b0 )

    xfm = np.dot( tkr_to_b0, reg )
    _, origin_ext = os.path.splitext( origin )

    if not origin_ext == 'vtk':
        # Call mris_convert
        new_origin = "%s_tmp.vtk" % output_prefix
        sp.check_call( 'mris_convert %s %s' % (origin, new_origin), shell=True )
        origin = new_origin
        del_vtk==True

    reader = tvtk.PolyDataReader( file_name=origin )
    vtk = reader.output
    reader.update()

    flip = np.identity(3)
    flip[0,0] *= -1.0;
    flip[1,1] *= -1.0;

    M = np.dot( flip, xfm[0:3,0:3] )
    O = np.dot( flip, xfm[0:3,3] )

    print M
    print O

    if use_paraview:
        phy_size = np.dot( ac[0:3,0:3], np.array(nib.load(target).get_data().shape) - [1,1,1] )
        #print phy_size
        #M[0,0:3]*=-1
        O[0]-=phy_size[0]*0.5
    
    for i,point in enumerate(vtk.points):
        vtk.points[i] = np.dot( M, point ) + O
    
    writer = tvtk.XMLPolyDataWriter( file_name='%s.vtp' % output_prefix, input=vtk )
#    writer = tvtk.PolyDataWriter( file_name='%s.vtk' % output_prefix, input=vtk )
    writer.update()

    savevtk( vtk, '%s.vtk' % output_prefix )

    if del_vtk:
        os.remove( origin )


def main():
    import sys
    parser = argparse.ArgumentParser( description="Transform freesurfer's surface to a targ space using an affine transformation" )
    parser.add_argument( '--targ', dest='target', type=lambda x: is_valid_file(parser,x), help='target volume (e.g. b=0 volume in dMRI)', required=True )
    parser.add_argument( '--orig', dest='origin', type=lambda x: is_valid_file(parser,x), help='original surface (e.g. surf/lh.pial outcome from freesurfer)', required=True  )
    parser.add_argument( '--o'   , dest='output_prefix', type=str, help='output prefix for filename (will be an XML Vtk file)' )
    parser.add_argument( '--tfm' , dest='transform', type=lambda x: is_valid_file(parser,x), help='transform file (e.g. register.dat from bbregister)', required=True  )
    parser.add_argument( '--paraview', dest='use_paraview', action='store_true', default=False, help='True if orientation should be fixed for Paraview rendering' )
    args = parser.parse_args( )

    if args.output_prefix == None:
        basename, ext = os.path.splitext( args.target )
        args.output_prefix = [ basename ]
        if ext == '.gz':
            args.output_prefix,_ = os.path.splitext( basename )
        args.output_prefix+= '_%s' % os.path.basename( args.origin )
    
    transform( args.origin, args.target, args.transform, args.output_prefix, args.use_paraview )

    return 0
    

def savevtk( vtk, fname, origin=(0.0,0.0,0.0) ):
    ncells = vtk.polys.number_of_cells
    npoints = len( vtk.points )
    nvalues = len( vtk.polys.data )
    cells = np.reshape( vtk.polys.data, ( ncells, int(nvalues/ncells) ) )
    now = dt.datetime.now()
    points = vtk.points + np.array(origin)

    if os.path.exists( fname ):
        os.remove( fname )

    try:
        with open(fname, 'a+') as f:
            f.write( '# vtk DataFile Version 2.0\n' ) # Write identification
            f.write( '%s, generated on %s\n' % (os.path.basename(fname),now.strftime("%Y-%m-%d %H:%M"))) # Write comment (up to 256 chars)
            f.write( 'ASCII\n' )
            f.write( 'DATASET POLYDATA\n' )
            f.write( 'POINTS %d float\n' % npoints ) 
            for row in points:
                f.write( '%.6f %.6f %.6f\n' % tuple(row) )
    
            f.write( 'POLYGONS %d %d\n' % ( ncells, nvalues ) )
            for row in cells:
                f.write( '3 %d %d %d\n' % tuple( row[1:] ) )
    except Exception as err:
        print err
        return 0

    return 1
    


if __name__ == '__main__':
    sys.exit(main())


