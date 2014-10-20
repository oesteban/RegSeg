# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# <h1>DWI Model Generation</h1>
# A model for testing the ACWE-Reg method is implemented here
# <ul>
#     <li>Author: phd@oscaresteban.es (Oscar Esteban)</li>
#     <li>Date: 2012/10/10</li>
#     <li>Version: 1.0</li>
# </ul>
# <h2>License</h2>
#
# <h2>Model Parameters</h2>
# Model parameters are estimated from a manual sample

# <codecell>

import os
model_path = os.path.abspath('Model1')

# <codecell>

import csv
import os
import numpy as np
from scipy.cluster.vq import *
from sklearn.covariance import EmpiricalCovariance, MinCovDet

col_l1 = 3

with open('./dti/dti_maps/20121011-ManualSampleOscar.csv', 'rU') as csvfile:
    dataReader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONE)
    dataReader.next()  # skip first row
    data = [[int(row[0]), float(row[col_l1]), float(
        row[col_l1 + 1]), float(row[col_l1 + 2])] for row in dataReader]
    csvfile.close()

getClass = lambda c: np.array(
    [[row[1], row[2], row[3]] for row in data if(row[0] == c)])

csf_data = getClass(1)
wm_data = getClass(2)
gm_data = getClass(3)

# <codecell>

for eval_data in (wm_data, gm_data, csf_data):
    robust_cov = MinCovDet().fit(eval_data)
    print robust_cov.location_
    print robust_cov.covariance_

# <codecell>

getTuple = lambda x: np.array([(row[0], row[1], row[2]) for row in x])

wm_c, wm_l = kmeans2(getTuple(wm_data), 3)
gm_c, gm_l = kmeans2(getTuple(gm_data), 3)
cs_c, cs_l = kmeans2(getTuple(csf_data), 3)

print wm_c
print gm_c
print cs_c

# <headingcell level=2>

# Tissue Fraction Distribution

# <codecell>

import nibabel as nib
import numpy as np
from scipy import ndimage
from pylab import *

# <codecell>

datashape = (101, 101, 101)
x, y, z = np.mgrid[0:datashape[0], 0:datashape[1], 0:datashape[2]]

# <codecell>


def ball(volsize, radius):
    assert radius < (volsize / 2)
    result = np.zeros(shape=(volsize, volsize, volsize), dtype=int)
    thres = radius ** 2
    center = np.floor(volsize / 2)
    for x in range(0, volsize + 1):
        for y in range(0, volsize + 1):
            for z in range(0, volsize + 1):
                val = (
                    (x - center) ** 2 + (y - center) ** 2 + (z - center) ** 2)
                if val < thres:
                    result[x, y, z] = 1
    return result

# <codecell>

# from
# http://stackoverflow.com/questions/9689173/shape-recognition-with-numpy-scipy-perhaps-watershed

import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt


def draw_circle(grid, x0, y0, radius):
    ny, nx = grid.shape
    y, x = np.ogrid[:ny, :nx]
    dist = np.hypot(x - x0, y - y0)
    grid[dist < radius] = True
    return grid

# <codecell>

data_wm = ball(101, 33)

# <codecell>

strel = ball(15, 6.5)

# <codecell>

data_gm = data_wm.copy()
for i in range(0, 1):
    data_gm = ndimage.binary_dilation(data_gm, structure=strel)
data_gm -= data_wm

# <codecell>

data_csf = np.zeros(shape=data_wm.shape)
data_csf[53, 57, 59] = 1
data_csf[45, 44, 40] = 1
for i in range(0, 2):
    data_csf = ndimage.binary_dilation(data_csf, structure=strel)
data_filled = np.copy(data_wm)
data_wm -= data_csf

# <codecell>

data_final = np.concatenate((data_wm[..., np.newaxis],
                             data_gm[..., np.newaxis],
                             data_csf[..., np.newaxis]), axis=3)
hdr = nib.spatialimages.Header(data_dtype=np.uint8, shape=data_final.shape)
niiHS = nib.Nifti1Image(data_final.astype(np.uint8), np.identity(4), hdr)
nib.save(niiHS, os.path.join(model_path, 'model.nii.gz'))

# <markdowncell>

# <h2>Contour extraction</h2>

# <codecell>


from tvtk.api import tvtk


def GenerateSurface(data, fname):
    grid = tvtk.ImageData(spacing=(1, 1, 1), origin=(0, 0, 0))
    grid.point_data.scalars = data.T.ravel()  # It wants fortran order???
    grid.point_data.scalars.name = 'scalars'
    grid.dimensions = data.shape
    iso = tvtk.ImageMarchingCubes(input=grid)
    w = tvtk.PolyDataWriter(
        input=iso.output, file_name=os.path.join(model_path, fname))
    w.write()


GenerateSurface(data_csf.astype(np.float), 'csf_prior.vtk')
GenerateSurface(data_filled.astype(np.float), 'wm_prior.vtk')
GenerateSurface(
    data_filled.astype(np.float) + data_gm.astype(np.float), 'gm_prior.vtk')
