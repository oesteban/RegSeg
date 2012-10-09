# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <markdowncell>

# <h1>DWI Model Generation</h1>
# A model for testing the ACWE-Reg method is implemented here
# Author: phd@oscaresteban.es (Oscar Esteban)
# Date: 2012/10/10
# Version: 1.0
# 
# License
# 
# <h2>Code</h2>

# <codecell>

import nibabel as nib
import numpy as np
from scipy import ndimage
from pylab import *

# <codecell>

datashape = ( 101, 101, 101 )
data = np.zeros( shape=datashape, dtype=int )

# <codecell>

data[ 50, 50, 50 ] = 1

# <codecell>

strel = np.zeros( shape=( 11, 11, 11 ), dtype=int )
strel[ 2:9,2:9, 2:9 ] = 1
strel[ 1:10,3:8, 3:8 ] = 1
strel[ 3:8,1:10, 3:8 ] = 1
strel[ 3:8, 3:8,1:10 ] = 1
strel[ :, :, 4 ]

# <codecell>

data2 = ndimage.binary_dilation( data, structure=strel )
data2 = ndimage.binary_dilation( data2, structure=strel )
data2 = ndimage.binary_dilation( data2, structure=strel )
data2 = ndimage.binary_dilation( data2, structure=strel )
im = imshow( data2[ :, :, 50 ] )

# <codecell>

data3 = ndimage.binary_dilation( data2, structure=strel )
data3 = ndimage.binary_dilation( data2, structure=strel )
im = imshow( data3[ :, :, 50 ] )

# <codecell>

data5 = np.zeros(shape=data.shape)
data5[ 53, 57, 59 ] = 1
data5[ 45, 44, 40 ] = 1
data5 = ndimage.binary_dilation( data5, structure=strel )

# <codecell>

data4 = np.zeros(shape=data.shape)
data4[ data3==1 ] = 3
data4[ data2==1 ] = 2
data4[ data5==1 ] = 1
im = imshow( data4[:,:, 40] )

# <codecell>

data_final = []
for i in range(1,4):
    data_class = np.zeros(shape=(101,101,101), dtype=np.uint8)
    data_class[data4 == i] = 1
    data_final.append( data_class )
#    figure()
#    imshow( data_class[:,:,50] )
np.shape( data_final )
data_final2 = np.concatenate((data_final[0][...,np.newaxis],data_final[1][...,np.newaxis],data_final[2][...,np.newaxis]),axis=3 )
#data_final2 = np.reshape(data_final, (101,101,101,3) )
#data_final2 = np.concatenate( ( [data_final[0]], [data_final[1]], [data_final[2]] ) )
np.shape( data_final2 )

# <codecell>

data_test = data_final2[:,:,:,1]
np.shape(data_test)
im=imshow( data_test[:,:,50] )

# <codecell>

niiHS = nib.Nifti1Image( data_final2, None, None )
nib.save( niiHS, 'model.nii.gz' )

# <codecell>


