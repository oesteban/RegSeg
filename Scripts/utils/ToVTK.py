import nibabel as nib
import numpy as np
from tvtk.api import tvtk
import os

from evtk.hl import imageToVTK 
# Dimensions 


fa = nib.load( os.path.join( os.path.abspath('Model1') , 'dtifit__FA.nii.gz') ).get_data()
(nx, ny, nz) = fa.shape 
ncells = nx * ny * nz 
npoints = (nx + 1) * (ny + 1) * (nz + 1) 

# Variables 
fa_data = np.array( fa, order = 'C') 

imageToVTK("./image", pointData = {"fa" : fa_data} )




#grid = tvtk.ImageData(spacing=(1,1,1), origin=(0,0,0))
#grid.point_data.scalars = fa.T.ravel().astype( np.float32 ) # It wants fortran order???
#grid.point_data.scalars.name = 'scalars'
#grid.dimensions = fa.shape
#w = tvtk.DataSetWriter(input=grid, file_name=os.path.join( os.path.abspath('Model1') , 'dtifit__FA.vtk' ) )
#w.write()
