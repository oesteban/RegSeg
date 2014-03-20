#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban
# @Date:   2014-03-11 15:28:16
# @Last Modified by:   oesteban
# @Last Modified time: 2014-03-20 12:41:31

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
    inputnode = pe.Node( niu.IdentityInterface(fields=['in_file','grid_size']), name='inputnode' )
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file', 'out_field', 'out_coeff' ]), name='outputnode' )
    distort = pe.Node( RandomBSplineDeformation(), name='bspline_field')

    wf.connect([
                 ( inputnode,  distort, [ ('in_file','in_file'), ('grid_size','grid_size')])
                ,( distort, outputnode, [ ('out_file','out_file'), ('out_field','out_field'), ('out_coeff','out_coeff')])
        ])

    return wf
