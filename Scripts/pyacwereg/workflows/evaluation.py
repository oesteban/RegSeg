#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: Oscar Esteban - code@oscaresteban.es
# @Date:   2014-03-12 16:59:14
# @Last Modified by:   Oscar Esteban
# @Last Modified time: 2014-03-12 17:39:09

import os
import os.path as op

import nipype.interfaces.io as nio              # Data i/o
import nipype.interfaces.utility as niu         # utility
import nipype.pipeline.engine as pe             # pypeline engine

from smri import prepare_smri
from distortion import bspline_deform

def bspline( name='BSplineEvaluation' ):
    """ A workflow to evaluate ACWEreg with random bspline deformations
    """
    wf = pe.Workflow( name=name )
    inputnode = pe.Node( niu.IdentityInterface( fields=[ 'subject_id', 'data_dir' ] ), name='inputnode' )
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file', 'out_field', 'out_coeff' ]), name='outputnode' )

    prep = prepare_smri()
    dist = bspline_deform()
    dist.inputs.inputnode.grid_size = 10

    wf.connect([
             ( inputnode, prep, [ ('subject_id','inputnode.subject_id'),('data_dir','inputnode.data_dir') ])
            ,( prep,      dist, [ ('outputnode.out_smri_brain','inputnode.in_file')])
            ,( dist, outputnode, [ ('outputnode.out_file','out_file'),('outputnode.out_field','out_field'),('outputnode.out_coeff','out_coeff')])
        ])

    return wf