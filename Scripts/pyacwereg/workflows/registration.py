#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban - code@oscaresteban.es
# @Date:   2014-03-28 20:38:30
# @Last Modified by:   oesteban
# @Last Modified time: 2014-04-04 12:41:35

import os
import os.path as op

import nipype.interfaces.io as nio              # Data i/o
import nipype.interfaces.utility as niu         # utility
import nipype.pipeline.engine as pe             # pipeline engine
import pyacwereg.nipype.interfaces as iface

def default_regseg( name='REGSEGDefault'):
    wf = pe.Workflow( name=name )

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_orig', 'in_dist', 'in_tpms', 'in_surf',
                        'in_mask' ]),
                        name='inputnode' )

    outputnode = pe.Node(niu.IdentityInterface(fields=['out_corr', 'out_tpms',
                         'out_surf', 'out_field', 'out_mask' ]),
                         name='outputnode' )

    # Registration
    regseg = pe.Node( iface.ACWEReg(), name="ACWERegistration" )
    regseg.inputs.iterations = [ 10, 10 ]
    regseg.inputs.descript_update = [ 5, 15 ]
    regseg.inputs.step_size = [ 0.5, 1.0 ]
    regseg.inputs.alpha = [ 0.0, 0.001 ]
    regseg.inputs.beta = [ 0.0, 0.01 ]
    regseg.inputs.grid_size = [ 6, 10 ]
    regseg.inputs.convergence_energy = [ True, True ]
    regseg.inputs.convergence_window = [ 5, 5 ]

    # Apply tfm to tpms
    applytfm = pe.Node( iface.FieldBasedWarp(), name="ApplyWarp" )

    # Connect
    wf.connect([
         ( inputnode,   regseg, [ ('in_surf','in_prior' ),
                                  ('in_dist','in_fixed' ) ])
        ,( inputnode, applytfm, [ ('in_tpms', 'in_file' ),
                                  ('in_mask', 'in_mask') ])
        ,( regseg,    applytfm, [ ('out_field','in_field')])
        ,( regseg,  outputnode, [ ('out_warped','out_corr'),
                                  ('out_field','out_field'),
                                  ('out_surfs', 'out_surf') ])
        ,( applytfm,outputnode, [ ('out_file', 'out_tpms'),
                                  ('out_mask', 'out_mask')])
    ])

    return wf

def identity_wf( name='Identity'):
    """ An identity workflow to check how ideal inverse transform
    affects final evaluation scores.
    """
    wf = pe.Workflow( name=name )

    def _invert_field( in_file, out_file=None ):
        import numpy as np
        import nibabel as nb
        import os.path as op

        im = nb.load( in_file )
        data = -1.0 * im.get_data()
        nii = nb.Nifti1Image( data, im.get_affine(), im.get_header() )
        if out_file is None:
            fname, ext = op.splitext( op.basename(in_file) )
            if ext == '.gz':
                fname, ext2 = op.splitext( fname )
                ext = ext2 + ext
            out_file = op.abspath( fname + '_inv' + ext )

        nb.save( nii, out_file )
        return out_file


    inputnode = pe.Node(niu.IdentityInterface(fields=['in_orig', 'in_dist', 'in_tpms', 'in_surf',
                        'in_mask', 'in_field' ]),
                        name='inputnode' )

    outputnode = pe.Node(niu.IdentityInterface(fields=['out_corr', 'out_tpms',
                         'out_surf', 'out_field', 'out_mask' ]),
                         name='outputnode' )

    # Invert field
    inv = pe.Node( niu.Function( input_names=['in_file'],
                   output_names=['out_file'], function=_invert_field ),
                   name='InvertField' )

    # Compute corrected images
    merge = pe.Node( niu.Merge(2), name='Merge' )
    split = pe.Node( niu.Split(splits=[2,3]), name='Split')

    # Apply tfm to tpms
    applytfm = pe.Node( iface.FieldBasedWarp(), name="ApplyWarp" )

    # Connect
    wf.connect([
         ( inputnode,       inv, [ ('in_field','in_file' ) ])
        ,( inputnode,     merge, [ ('in_dist', 'in1'), ('in_tpms', 'in2' ) ])
        ,( inputnode,  applytfm, [ ('in_mask', 'in_mask'),
                                   ('in_surf', 'in_surf') ])
        ,( merge,      applytfm, [ ('out', 'in_file' )])
        ,( inv,        applytfm, [ ('out_file', 'in_field')])
        ,( applytfm,      split, [ ('out_file', 'inlist')])
        ,( split,    outputnode, [ ('out1', 'out_corr' ),
                                   ('out2', 'out_tpms' )])
        ,( inv,      outputnode, [ ('out_file','out_field') ])
        ,( applytfm, outputnode, [ ('out_surf', 'out_surf'),
                                   ('out_mask', 'out_mask')])
    ])

    return wf
