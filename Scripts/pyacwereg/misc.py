#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
misc.py - Miscelaneous helpers for ACWEReg

Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban),
                    with Biomedical Image Technology, UPM (BIT-UPM)
All rights reserved.
This file is part of ACWEReg.

"""

__author__ = "Oscar Esteban"
__copyright__ = "Copyright 2013, Biomedical Image Technologies (BIT), \
                 Universidad Politécnica de Madrid"
__credits__ = ["Oscar Esteban"]
__license__ = "FreeBSD"
__version__ = "0.1"
__maintainer__ = "Oscar Esteban"
__email__ = "code@oscaresteban.es"
__status__ = "Prototype"


import os.path as op
import numpy as np
import nibabel as nb
from math import sqrt
from scipy import ndimage


def sort_surfs(surfs):
    import os.path as op

    if not isinstance(surfs, list):
        return surfs
    if len(surfs) == 1:
        return surfs

    out = [None] * len(surfs)
    inames = [s.lower() for s in op.basename(surfs)]

    for fname, iname in zip(surfs, surfs):
        index = -1
        if 'white' in iname:
            if ('lh' in iname) or ('.l.' in iname):
                index = 0
            elif ('rh' in iname) or ('.r' in iname):
                index = 1
        elif 'pial' in iname:
            if ('lh' in iname) or ('.l.' in iname):
                index = 2
            elif ('rh' in iname) or ('.r' in iname):
                index = 3
        out.insert(index, fname)

    return out


def ball(volsize, radius, dims=3):
    volsize = np.array(volsize)
    if volsize.ndim == 0:
        volsize = np.ones(shape=(dims)) * volsize

    assert np.all(volsize > 0)
    assert float(radius) < (float(volsize.min()) * 0.5)

    result = np.zeros(shape=tuple(volsize), dtype=int)
    center = (volsize - 1) * 0.5

    for x in range(0, int(volsize[0])):
        for y in range(0, int(volsize[1])):
            for z in range(0, int(volsize[2])):
                if np.linalg.norm(center - [x, y, z]) < radius:
                    result[x, y, z] = 1
    return result


def gen_noise(image, mask=None, snr_db=10.0):
    """
    Generates a copy of an image with a certain amount of
    added gaussian noise (rayleigh for background in mask)
    """
    snr = np.power(10.0, snr_db / 10.0)
    noise = np.random.normal(size=image.shape)
    bg_noise = np.random.rayleigh(size=image.shape)

    if mask is None:
        mask = np.ones_like(image)

    im_scaled = noise.std() / image.std() * (sqrt(snr)) * image
    im_noise = np.zeros_like(image)

    im_noise[mask > 0] = im_scaled[mask > 0] + noise[mask > 0]
    im_noise[mask == 0] = bg_noise[mask == 0]

    return im_noise


def genNiftiVol(data, dtype=np.uint8):
    shape = np.array(np.shape(data), dtype=np.float32)
    if np.ndim(data) > 3:
        shape = shape[1:]

    affine = np.identity(4)
    affine[0:3, 3] = -0.5 * shape[:3]
    hdr = nb.Nifti1Header()
    hdr.set_data_dtype(np.uint8)
    hdr.set_xyzt_units('mm')
    hdr.set_data_dtype(dtype)

    hdr['data_type'] = 2
    hdr['qform_code'] = 2  # aligned
    hdr['sform_code'] = 1  # scanner
    hdr['scl_slope'] = 1.0
    hdr['regular'] = np.array('r', dtype='|S1')
    pixdim = np.ones(shape=(8,))
    pixdim[4:] = 0
    hdr['pixdim'] = pixdim

    if np.ndim(data) > 3:
        hdr.set_xyzt_units('mm', 'sec')
        # hdr['xyzt_units'] = 2 + 8  # mm + sec
        pixdim[4] = 1
        hdr['pixdim'] = pixdim
        nii_array = []
        for im in data:
            nii_array.append(nb.Nifti1Image(im.astype(dtype), affine, hdr))
        nii = nb.concat_images(nii_array)

    else:
        nii = nb.Nifti1Image(data.astype(dtype), affine, hdr)
    return nii


def genBall(datashape=(101, 101, 101), radius=17, cortex=True):
    wm = ball(datashape, radius)

    if cortex:
        ball2 = ball(11, 4.4)
        gm = ndimage.binary_dilation(wm, structure=ball2).astype(np.uint8) - wm
        bg = np.ones_like(wm) - (gm + wm)
        return [bg, wm, gm]
    else:
        bg = np.ones_like(wm) - wm
        return [bg, wm]


def genGyrus(datashape=(101, 101, 101), radius=35, cortex=True):
    modelbase = ball(datashape, radius)
    center_pix = ((np.array(datashape) - 1) * 0.5).astype(np.uint8)
    displ_pix = ((np.array(datashape) - 1) * 0.25).astype(np.uint8)
    modelbase[center_pix[0], center_pix[1]:, :] = 0
    modelbase[:displ_pix[0], center_pix[1], :] = 0
    ball1 = ball(11, 4.5)
    wm = ndimage.binary_opening(ndimage.binary_erosion(
        modelbase, structure=ball1).astype(np.uint8),
        structure=ball1).astype(np.uint8)

    if cortex:
        ball2 = ball(11, 4.4)
        gm = ndimage.binary_dilation(wm, structure=ball2).astype(np.uint8) - wm
        bg = np.ones_like(modelbase) - (gm + wm)
        return [bg, wm, gm]
    else:
        bg = np.ones_like(modelbase) - wm
        return [bg, wm]


def genBox(datashape=(101, 101, 101), coverage=0.4, cortex=True):
    modelbase = np.zeros(shape=datashape)
    extent = np.around(coverage * np.array(datashape))
    padding = np.around(0.5 * (np.array(datashape) - extent))
    end = np.array(datashape) - padding
    modelbase[padding[0]:end[0], padding[1]:end[1], padding[2]:end[2]] = 1

    ball1 = ball(11, 4.5)
    wm = ndimage.binary_opening(ndimage.binary_erosion(
        modelbase, structure=ball1).astype(np.uint8),
        structure=ball1).astype(np.uint8)

    if cortex:
        ball2 = ball(11, 4.4)
        gm = ndimage.binary_dilation(wm, structure=ball2).astype(np.uint8) - wm
        bg = np.ones_like(modelbase) - (gm + wm)
        return [bg, wm, gm]
    else:
        bg = np.ones_like(modelbase) - wm
        return [bg, wm]


def genL(datashape=(101, 101, 101), cortex=True):
    modelbase = np.zeros(shape=datashape)
    center = np.around(0.5 * np.array(datashape))
    extent = np.around(0.4 * np.array(datashape))
    padding = np.around(0.5 * (np.array(datashape) - extent))
    end = np.array(datashape) - padding
    modelbase[padding[0]:end[0], padding[1]:end[1], padding[2]:end[2]] = 1
    modelbase[center[0]:end[0], center[1]:end[1], center[2]:end[2]] = 0

    ball1 = ball(11, 4.5)
    wm = ndimage.binary_opening(ndimage.binary_erosion(
        modelbase, structure=ball1).astype(np.uint8), structure=ball1).astype(np.uint8)

    if cortex:
        ball2 = ball(11, 4.4)
        gm = ndimage.binary_dilation(wm, structure=ball2).astype(np.uint8) - wm
        bg = np.ones_like(modelbase) - (gm + wm)
        return [bg, wm, gm]
    else:
        bg = np.ones_like(modelbase) - wm
        return [bg, wm]


def genShape(name, datashape=(101, 101, 101), cortex=True):
    if name == 'box':
        return genBox(datashape=datashape, cortex=cortex)
    elif name == 'L':
        return genL(datashape=datashape, cortex=cortex)
    elif name == 'ball':
        return genBall(datashape=datashape, cortex=cortex)
    elif name == 'gyrus':
        return genGyrus(datashape=datashape, cortex=cortex)
    else:
        return genBox()


def genContrast(model, values):
    assert(len(model) > 1)
    assert((len(model) - 1) <= len(values))

    if((len(model) - 1) < len(values)):
        values = values[0:len(model)]

    contrast = np.zeros_like(model[0])
    for c, v in zip(model[1:], values):
        contrast = contrast + c * v
    return contrast


def draw_circle(grid, x0, y0, radius):
    """
    http://stackoverflow.com/questions/9689173/shape-recognition-with-numpy-scipy-perhaps
    """
    ny, nx = grid.shape
    y, x = np.ogrid[:ny, :nx]
    dist = np.hypot(x - x0, y - y0)
    grid[dist < radius] = True
    return grid


def genSurface(data, fname):
    from tvtk.api import tvtk
    grid = tvtk.ImageData(spacing=(1, 1, 1), origin=(0, 0, 0))
    grid.point_data.scalars = data.T.ravel()  # It wants fortran order???
    grid.point_data.scalars.name = 'scalars'
    grid.dimensions = data.shape
    iso = tvtk.ImageMarchingCubes(input=grid)
    w = tvtk.PolyDataWriter(
        input=iso.output, file_name=os.path.join(model_path, fname))
    w.write()


def nii2vtk(in_file, out_file=None):
    from evtk.hl import imageToVTK
    nii = nb.load(in_file)
    data = np.array(nii.get_data(), order='C')

    if out_file is None:
        out_file, ext = op.splitext(op.basename(in_file))
        if ext == '.gz':
            out_file, _ = op.splitext(out_file)

    out_file = op.abspath(out_file)
    imageToVTK(out_file, pointData={'scalar': data})
    return out_file
