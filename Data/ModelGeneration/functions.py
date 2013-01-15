import os
import numpy as np
from tvtk.api import tvtk
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

def GenerateSurface( data, fname, use_smooth=True, spacing=(1.0,1.0,1.0), origin=(0.0,0.0,0.0) ):
    grid = tvtk.ImageData(spacing=spacing, origin=origin, dimensions=data.shape)
    grid.point_data.scalars = data.ravel().astype( np.uint8 )
    grid.point_data.scalars.name = 'scalars'
    iso = tvtk.ImageMarchingCubes(input=grid)
    iso.update()
    contour = iso.output
    if use_smooth:
        vtkmeshclean = tvtk.CleanPolyData( input=iso.output )
        vtkmeshclean.update()
        smooth = tvtk.SmoothPolyDataFilter( input=vtkmeshclean.output, number_of_iterations=20 )
        smooth.update()
        vtkmeshclean2 = tvtk.CleanPolyData( input=smooth.output )
        vtkmeshclean2.update()
        contour = vtkmeshclean2.output
    contour.point_data.scalars = np.ones( np.shape( contour.points )[0] )
    contour.point_data.scalars.name = 'scalars'
    w = tvtk.PolyDataWriter(input=contour, file_name=fname )
    w.write()
    return contour


def FilterMesh(surf, sl, axis=2 ):
    points = np.array(surf.points, dtype=np.float32 )
    filtered = points[(points[:,axis]>(sl-0.25)) & (points[:,axis]<(sl+0.25)) ]
    filtered[:,axis] = sl
    uniq = np.array(list(set(tuple(p) for p in filtered)))
    mesh = tvtk.PolyData( points=uniq )
    mesh.point_data.scalars = np.zeros( np.shape(uniq)[0], 'f' )
    mesh.point_data.scalars.name = 'scalars'
#    surf = tvtk.SurfaceReconstructionFilter(input=mesh)
#    surf.update()
#    cont = tvtk.ContourFilter( input=mesh )
#    cont.update()
#    rev = tvtk.ReverseSense( input=cont.output, reverse_cells=True, reverse_normals=True)
#    rev.update()
#    mesh = cont.output
    return mesh

def ComputeContourFunction( image, surf, sl, axis=2, fig=None, color='r' ):
    dim = ( image.get_shape()[0],image.get_shape()[1],image.get_shape()[2] )
    dim2d = tuple(np.delete( np.array(dim), axis ))
    mesh = FilterMesh( surf, sl, axis )
    ps = np.array( mesh.points, dtype=np.float32 )[:,0:2]
    d = np.ones( ps.shape[0] ).astype( np.float32 )    
    xi,yi = np.mgrid[0:dim2d[0],0:dim2d[1]] 
    interpolated = griddata( ps, d, (xi,yi), method='linear', fill_value=0.0 )
    return interpolated
