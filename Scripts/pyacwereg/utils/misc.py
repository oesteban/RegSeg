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


import numpy as np
import nibabel as nib
from math import sqrt


def ball( volsize, radius, dims=3):
    volsize = np.array( volsize )
    if volsize.ndim == 0:
        volsize = np.ones( shape=(dims) ) * volsize

    assert np.all( volsize>0 )
    assert float(radius) < ( float( volsize.min() ) * 0.5 )

    result = np.zeros( shape=tuple(volsize), dtype=int )
    center = (volsize-1) * 0.5

    for x in range(0,int(volsize[0])):
        for y in range(0,int(volsize[1])):
            for z in range(0,int(volsize[2])):
                if  np.linalg.norm( center - [x,y,z] ) < radius:
                    result[x,y,z] = 1
    return result


def gen_noise( image, mask=None, snr_db=10.0 ):
    """
    Generates a copy of an image with a certain amount of
    added gaussian noise (rayleigh for background in mask)
    """
    snr = np.power( 10.0, snr_db/10.0 )
    noise = np.random.normal( size=image.shape )
    bg_noise = np.random.rayleigh( size=image.shape )

    if mask is None:
        mask = np.ones_like( image )

    im_scaled = noise.std() / image.std() *(sqrt(snr))* image;
    im_noise = np.zeros_like( image )

    im_noise[mask>0] = im_scaled[mask>0] + noise[mask>0]
    im_noise[mask==0] = bg_noise[mask==0]

    return im_noise


def normalize_tpms( in_files, in_mask=None, out_files=[] ):
    import nibabel as nib
    import numpy as np
    import os.path as op

    in_files = np.atleast_1d( in_files ).tolist()

    if len( out_files )!=len(in_files):
        for i,finname in enumerate( in_files ):
            fname,fext = op.splitext( op.basename( finname ) )
            if fext == '.gz':
                fname,fext2 = op.splitext( fname )
                fext = fext2 + fext

            out_file = op.abspath( fname+'_norm'+fext )
            out_files+= [ out_file ]

    imgs = [ nib.load(fim) for fim in in_files ]

    if len(in_files)==1:
        img_data = imgs[0].get_data()
        img_data[img_data>0.0] = 1.0
        hdr = imgs[0].get_header().copy()
        hdr['data_type']= 16
        hdr.set_data_dtype( 'float32' )
        nib.save( nib.Nifti1Image( img_data.astype(np.float32), imgs[0].get_affine(), hdr ), out_files[0] )
        return out_files[0]

    img_data = np.array( [ im.get_data() for im in imgs ] ).astype( 'f32' )
    #img_data[img_data>1.0] = 1.0
    img_data[img_data<0.0] = 0.0
    weights = np.sum( img_data, axis=0 )

    msk = np.ones_like( imgs[0].get_data() )
    msk[ weights<= 0 ] = 0

    if not in_mask is None:
        msk = nib.load( in_mask ).get_data()
        msk[ msk<=0 ] = 0
        msk[ msk>0 ] = 1

    msk = np.ma.masked_equal( msk, 0 )


    for i,out_file in enumerate( out_files ):
        data = np.ma.masked_equal( img_data[i], 0 )
        probmap = data / weights
        hdr = imgs[i].get_header().copy()
        hdr['data_type']= 16
        hdr.set_data_dtype( 'float32' )
        nib.save( nib.Nifti1Image( probmap.astype(np.float32), imgs[i].get_affine(), hdr ), out_file )

    return out_files


def genNiftiVol( data, dtype=np.uint8 ):
    import numpy as np
    import nibabel as nb

    shape = np.array(np.shape( data ), dtype=np.float32 )
    if np.ndim( data ) > 3:
        shape = shape[1:]

    affine = np.identity( 4 )
    affine[0:3,3] = -0.5 * shape
    hdr = nb.Nifti1Header()
    hdr.set_data_dtype( np.uint8 )
    hdr['xyzt_units'] = 2 # mm.

    hdr['data_type'] = 2
    hdr['qform_code'] = 2 # aligned
    hdr['sform_code'] = 1 # scanner
    hdr['scl_slope'] = 1.0
    hdr['regular'] = np.array('r', dtype='|S1')
    pixdim = np.ones( shape=(8,) )
    pixdim[4:] = 0
    hdr['pixdim'] = pixdim

    if np.ndim( data ) > 3:
        hdr['xyzt_units'] = 2 + 8 # mm + sec
        pixdim[4] = 1
        hdr['pixdim'] = pixdim
        nii_array = []
        for im in data:
            nii_array.append( nb.Nifti1Image( im.astype(dtype), affine, hdr ) )
        nii = nb.concat_images( nii_array )

    else:
        nii = nb.Nifti1Image( data.astype( dtype ), affine, hdr )
    return nii

def genBall(datashape=( 101,101,101 ), radius=17, cortex=True):
    import pyacwereg.utils.misc as misc
    import scipy.ndimage as ndimage
    import numpy as np

    wm = ball( datashape, radius )

    if cortex:
        ball2 = ball(11,4.4)
        gm = ndimage.binary_dilation( wm, structure=ball2 ).astype( np.uint8 ) - wm
        bg = np.ones_like( wm ) - (gm + wm)
        return [ bg, wm, gm ]
    else:
        bg = np.ones_like( wm ) - wm
        return [ bg, wm ]

def genGyrus(datashape=(101,101,101), radius=35, cortex=True):
    import pyacwereg.utils.misc as misc
    import scipy.ndimage as ndimage
    import numpy as np

    modelbase = ball( datashape, radius )
    center_pix = ((np.array( datashape )-1)*0.5).astype(np.uint8)
    modelbase[center_pix[0],:, center_pix[2]: ] = 0
    ball1 = ball(11,4.5)
    wm = ndimage.binary_opening( ndimage.binary_erosion( modelbase, structure=ball1 ).astype( np.uint8 ), structure=ball1 ).astype( np.uint8 )

    if cortex:
        ball2 = ball(11,4.4)
        gm = ndimage.binary_dilation( wm, structure=ball2 ).astype( np.uint8 ) - wm
        bg = np.ones_like( modelbase ) - (gm + wm)
        return [ bg, wm, gm ]
    else:
        bg = np.ones_like( modelbase ) - wm
        return [ bg, wm ]

def genBox( datashape=(101,101,101), coverage=0.4, cortex=True ):
    import pyacwereg.utils.misc as misc
    import scipy.ndimage as ndimage
    import numpy as np

    modelbase = np.zeros( shape=datashape )
    extent = np.around(  coverage * np.array( datashape ) )
    padding = np.around( 0.5 * (np.array( datashape ) - extent) )
    end = np.array( datashape ) - padding
    modelbase[padding[0]:end[0],padding[1]:end[1],padding[2]:end[2]] = 1

    ball1 = ball(11,4.5)
    wm = ndimage.binary_opening( ndimage.binary_erosion( modelbase, structure=ball1 ).astype( np.uint8 ), structure=ball1 ).astype( np.uint8 )

    if cortex:
        ball2 = ball(11,4.4)
        gm = ndimage.binary_dilation( wm, structure=ball2 ).astype( np.uint8 ) - wm
        bg = np.ones_like( modelbase ) - (gm + wm)
        return [ bg, wm, gm ]
    else:
        bg = np.ones_like( modelbase ) - wm
        return [ bg, wm ]

def genL( datashape=(101,101,101), cortex=True ):
    import scipy.ndimage as ndimage
    import numpy as np

    modelbase = np.zeros( shape=datashape )
    center = np.around(  0.5 * np.array( datashape ) )
    extent = np.around(  0.4 * np.array( datashape ) )
    padding = np.around( 0.5 * (np.array( datashape ) - extent) )
    end = np.array( datashape ) - padding
    modelbase[padding[0]:end[0],padding[1]:end[1],padding[2]:end[2]] = 1
    modelbase[center[0]:end[0],center[1]:end[1],center[2]:end[2]] = 0

    ball1 = ball(11,4.5)
    wm = ndimage.binary_opening( ndimage.binary_erosion( modelbase, structure=ball1 ).astype( np.uint8 ), structure=ball1 ).astype( np.uint8 )

    if cortex:
        ball2 = ball(11,4.4)
        gm = ndimage.binary_dilation( wm, structure=ball2 ).astype( np.uint8 ) - wm
        bg = np.ones_like( modelbase ) - (gm + wm)
        return [ bg, wm, gm ]
    else:
        bg = np.ones_like( modelbase ) - wm
        return [ bg, wm ]

def genShape( name, cortex=True ):
    if name == 'box':
        return genBox( cortex=cortex )
    elif name == 'L':
        return genL( cortex=cortex )
    elif name == 'ball':
        return genBall( cortex=cortex )
    elif name == 'gyrus':
        return genGyrus( cortex=cortex )
    else:
        return genBox()

def genContrast( model, values ):
    assert( len( model ) > 1 )
    assert( (len(model)-1) <= len(values) )

    if( (len(model)-1) < len(values) ):
        values = values[0:len(model)]


    contrast = np.zeros_like( model[0] )
    for c,v in zip( model[1:], values ):
        contrast = contrast + c * v
    return contrast
