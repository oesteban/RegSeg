#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban
# @Date:   2014-03-11 15:28:16
# @Last Modified by:   oesteban
# @Last Modified time: 2014-04-15 09:07:55

import os
import os.path as op

import nipype.interfaces.io as nio              # Data i/o
import nipype.interfaces.utility as niu         # utility
import nipype.pipeline.engine as pe             # pipeline engine
import nipype.interfaces.fsl as fsl             # fsl
import nipype.interfaces.freesurfer as fs       # freesurfer
import nipype.interfaces.ants as ants           # ANTS
import nipype.pipeline.engine as pe

from epi import epi_deform,fieldmap_preparation
from pyacwereg.nipype.interfaces import RandomBSplineDeformation
from pyacwereg.utils.misc import normalize_tpms as normalize

def isbi_workflow( name='ISBI2014' ):
    workflow = pe.Workflow(name=name)
    # Setup i/o
    inputnode = pe.Node( niu.IdentityInterface(fields=['subject_id','in_fmap_mag','in_fmap_pha','in_t1w', 'in_t2w', 'fs_subjects_dir','te_incr', 'echospacing','enc_dir' ]), name='inputnode' )
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file', 'out_vsm', 'out_mask', 'out_tpms' ]), name='outputnode' )

    # Set-up internal workflows
    prepare = fieldmap_preparation()
    distort = epi_deform()

    workflow.connect([
                         ( inputnode,  prepare, [ ('subject_id','inputnode.subject_id'),('in_fmap_mag','inputnode.in_fmap_mag'),
                                                  ('in_fmap_pha','inputnode.in_fmap_pha'),('in_t1w','inputnode.in_t1w'),
                                                  ('in_t2w','inputnode.in_t2w'),('fs_subjects_dir','inputnode.fs_subjects_dir') ])
                        ,( inputnode,  distort, [ ('te_incr','inputnode.te_incr'),('echospacing','inputnode.echospacing'),('enc_dir','inputnode.enc_dir') ])
                        ,( prepare,    distort, [ ('outputnode.out_smri','inputnode.in_file'),('outputnode.out_fmap_pha','inputnode.in_pha'),
                                                  ('outputnode.out_fmap_mag','inputnode.in_mag'),('outputnode.out_mask','inputnode.in_mask'),
                                                  ('outputnode.out_tpms','inputnode.in_tpms') ])
                        ,( distort, outputnode, [ ('outputnode.out_file','out_file'),('outputnode.out_vsm','out_vsm'),
                                                  ('outputnode.out_mask','out_mask'),('outputnode.out_tpms','out_tpms') ])
                    ])

    return workflow

def bspline_deform( name='BSplineDistortion' ):
    """ A nipype workflow to produce bspline-based deformation fields """
    wf = pe.Workflow( name=name )
    inputnode = pe.Node( niu.IdentityInterface(fields=['in_file','in_tpms','in_mask','in_surfs','grid_size']),
                         name='inputnode' )
    outputnode = pe.Node( niu.IdentityInterface(fields=['out_file','out_tpms', 'out_surfs',
                                                       'out_field', 'out_coeff','out_mask' ]),
                          name='outputnode' )

    merge = pe.Node( niu.Merge(2), name='MergeInputs')

    distort = pe.Node( RandomBSplineDeformation(), name='bspline_field')

    split = pe.Node( niu.Split(splits=[2,3]), name='SplitOutputs' )

    norm_tpms = pe.Node( niu.Function( input_names=['in_files','in_mask'], output_names=['out_files'], function=normalize ), name='Normalize' )

    wf.connect([
       ( inputnode,    distort, [ ('in_surfs','in_surfs'),
                                  ('grid_size','grid_size'),
                                  ('in_mask', 'in_mask') ])
      ,( inputnode,      merge, [ ('in_file', 'in1'), ('in_tpms','in2') ])
      ,( merge,        distort, [ ('out', 'in_file' )])
      ,( distort,        split, [ ('out_file', 'inlist') ])
      ,( split,      norm_tpms, [ ('out2','in_files' )])
      ,( distort,    norm_tpms, [ ('out_mask', 'in_mask' )])
      ,( split,     outputnode, [ ('out1', 'out_file') ])
      ,( norm_tpms, outputnode, [ ('out_files','out_tpms')])
      ,( distort,   outputnode, [ ('out_surfs', 'out_surfs'),
                                  ('out_field','out_field'),
                                  ('out_coeff','out_coeff'),
                                  ('out_mask','out_mask')])
        ])

    return wf
