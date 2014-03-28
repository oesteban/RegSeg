#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: Oscar Esteban - code@oscaresteban.es
# @Date:   2014-03-12 16:59:14
# @Last Modified by:   oesteban
# @Last Modified time: 2014-03-28 13:06:28

import os
import os.path as op

import nipype.interfaces.io as nio              # Data i/o
import nipype.interfaces.utility as niu         # utility
import nipype.pipeline.engine as pe             # pipeline engine
import pyacwereg.nipype.interfaces as iface

from smri import prepare_smri
from distortion import bspline_deform

def bspline( name='BSplineEvaluation' ):
    """ A workflow to evaluate ACWEreg with random bspline deformations
    """
    wf = pe.Workflow( name=name )
    inputnode = pe.Node( niu.IdentityInterface( fields=[ 'subject_id', 'data_dir','grid_size' ] ), name='inputnode' )
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file', 'out_tpms', 'out_surfs',
                                               'out_field', 'out_coeff' ]), name='outputnode' )

    prep = prepare_smri()
    dist = bspline_deform()
    regseg = pe.Node( iface.ACWEReg(), name="ACWERegistration" )
    regseg.inputs.iterations = [ 10, 10 ]
    regseg.inputs.descript_update = [ 5, 15 ]
    regseg.inputs.step_size = [ 0.5, 1.0 ]
    regseg.inputs.alpha = [ 0.0, 0.001 ]
    regseg.inputs.beta = [ 0.0, 0.01 ]
    regseg.inputs.grid_size = [ 6, 10 ]
    regseg.inputs.convergence_energy = [ True, True ]
    regseg.inputs.convergence_window = [ 5, 5 ]

    wf.connect([
             ( inputnode,  prep, [ ('subject_id','inputnode.subject_id'),('data_dir','inputnode.data_dir') ])
            ,( inputnode,  dist, [ ('grid_size', 'inputnode.grid_size')])
            ,( prep,       dist, [ ('outputnode.out_smri_brain','inputnode.in_file'),
                                   ('outputnode.out_surfs', 'inputnode.in_source')
                                   ('outputnode.out_tpms', 'inputnode.in_tpms')])
            ,( prep,     regseg, [ ('outputnode.out_surfs','in_prior' ) ])
            ,( dist,     regseg, [ ('outputnode.out_file','in_fixed' ) ])
            ,( dist, outputnode, [ ('outputnode.out_file','out_file'),('outputnode.out_field','out_field'),('outputnode.out_coeff','out_coeff')])
        ])

    return wf
