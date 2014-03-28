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

    imgs = [ nib.load(fim) for fim in in_files ]
    img_data = np.array( [ im.get_data() for im in imgs ] ).astype( 'f32' )
    img_data[img_data>1.0] = 1.0
    img_data[img_data<0.0] = 0.0
    weights = np.sum( img_data, axis=0 )

    img_data[0][weights==0] = 1.0
    weights[weights==0] = 1.0

    msk = np.ones_like( imgs[0].get_data() )

    if not in_mask is None:
        msk = nib.load( in_mask ).get_data()
        msk[ msk<=0 ] = 0
        msk[ msk>0 ] = 1


    if len( out_files )==0:
        for i,finname in enumerate( in_files ):
            fname,fext = op.splitext( op.basename( finname ) )
            if fext == '.gz':
                fname,fext2 = op.splitext( fname )
                fext = fext2 + fext

            out_file = op.abspath( fname+'_norm'+fext )
            out_files+= [ out_file ]


    for i,out_file in enumerate( out_files ):
            data = img_data[i] / weights
            data = data * msk
            hdr = imgs[i].get_header().copy()
            hdr['data_type']= 16
            hdr.set_data_dtype( 'float32' )
            nib.save( nib.Nifti1Image( data, imgs[i].get_affine(), hdr ), out_file )

    return out_files
