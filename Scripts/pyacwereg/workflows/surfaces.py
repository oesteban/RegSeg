#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
surface.py: Implements workflows relative to surfaces manipulation
            for ACWEReg

Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban)
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

import nipype.pipeline.engine as pe
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as niu
import nipype.interfaces.io as nio

import sys
sys.path.append('/home/oesteban/workspace/ACWE-Reg/Scripts/utils')
import pyacwereg.utils.fsb_transform as fsbt


def mris_convert( in_file, extension="vtk", out_file=None ):
    import subprocess as sp
    import numpy as np
    import os.path as op

    if out_file is None:
        out_file = op.abspath('./%s_conv' % op.basename( in_file ) )

    if not extension is None and extension!="":
        out_file = out_file + ".%s" % extension

    cmd = "mris_convert %s %s" % (in_file, out_file )
    proc = sp.Popen( cmd, stdout=sp.PIPE, shell=True )

    return out_file

def binary2contour_workflow( name="contour" ):
    pipeline = pe.Workflow( name=name )
    inputnode = pe.Node(niu.IdentityInterface( fields=['in_file','in_norm','out_folder' ] ), name='inputnode' )
    outputnode = pe.Node(niu.IdentityInterface( fields=['out_file'] ), name='outputnode' )

    #convertmgz = pe.Node( fs.MRIConvert( out_type='mgz' ), name='convertmgz' )
    pretess = pe.Node( niu.Function( input_names=['in_file', 'in_norm' ], output_names=['out_file'],
                                    function=fsbt.MRIPreTess ), name='pretess' )
    tess = pe.Node( fs.MRITessellate(label_value=1,use_real_RAS_coordinates=True), name='tess' )
    smooth = pe.Node( fs.SmoothTessellation(disable_estimates=True ), name='mris_smooth' )
    toVtk = pe.Node( niu.Function( input_names=['in_file'], output_names=['out_file'], function=mris_convert ), name="tovtk" )
    ds = pe.Node( nio.DataSink( container='', parameterization=False ), name='write' )


    pipeline.connect( [

                        (inputnode, pretess, [ ('in_file','in_file'), ('in_norm', 'in_norm' ) ] )
                       ,(pretess,tess,       [ ('out_file', 'in_file' ) ])
                       ,(tess,smooth,        [ ('surface', 'in_file' ) ])
                       ,(smooth, toVtk,      [ ('surface', 'in_file' ) ])
                       ,(toVtk, ds,          [ ('out_file', 'surfs' ) ])
                       ,(inputnode,ds,       [ ('out_folder','base_directory') ])
                       ,(ds, outputnode,     [ ('out_file','out_file') ])
                       ])

    return pipeline
