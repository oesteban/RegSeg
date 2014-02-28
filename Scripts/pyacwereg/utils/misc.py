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


