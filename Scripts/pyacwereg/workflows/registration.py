#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban - code@oscaresteban.es
# @Date:   2014-03-28 20:38:30
# @Last Modified by:   oesteban
# @Last Modified time: 2014-10-23 12:51:15

import os
import os.path as op

import nipype.pipeline.engine as pe             # pipeline engine
from nipype.interfaces import io as nio              # Data i/o
from nipype.interfaces import utility as niu         # utility

from pyacwereg.interfaces.acwereg import ACWEReg
from pyacwereg.interfaces.warps import FieldBasedWarp, InverseField


def default_regseg(name='REGSEGDefault'):
    wf = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
                        fields=['in_orig', 'in_dist', 'in_tpms', 'in_surf',
                                'in_mask', 'grid_size']), name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
                         fields=['out_corr', 'out_tpms',
                                 'out_surf', 'out_field', 'out_mask']),
                         name='outputnode')

    # Registration
    # Good config for box phantom (2014/04/21): [ -a 0.0 -b 0.0 -u 20 -g 6 -i
    # 500 -s 1.0]
    regseg = pe.Node(ACWEReg(), name="ACWERegistration")
    regseg.inputs.iterations = [500, 500, 500]
    #regseg.inputs.descript_update = [20]
    regseg.inputs.step_size = [.001, .1, .1]
    regseg.inputs.alpha = [0.0, 0.0, 0.0]
    regseg.inputs.beta = [0.0, 0.0, 0.0]
    regseg.inputs.grid_size = [4, 5, 6]
    regseg.inputs.convergence_energy = [True] * 3
    regseg.inputs.convergence_window = [50, 25, 15]
    regseg.inputs.f_smooth = [2.0, 1.0, None]

    # Apply tfm to tpms
    applytfm = pe.MapNode(FieldBasedWarp(), name="ApplyWarp",
                          iterfield=['in_file'])

    # Connect
    wf.connect([
        (inputnode,   regseg, [('in_surf', 'in_prior'),
                               ('in_dist', 'in_fixed')]),
        (inputnode, applytfm, [('in_tpms', 'in_file'),
                               ('in_mask', 'in_mask')]),
        (regseg,    applytfm, [('out_field', 'in_field')]),
        (regseg,  outputnode, [('out_warped', 'out_corr'),
                               ('out_field', 'out_field'),
                               ('out_surfs', 'out_surf')]),
        (applytfm, outputnode, [('out_file', 'out_tpms'),
                                ('out_mask', 'out_mask')])
    ])

    return wf


def identity_wf(name='Identity', n_tissues=3):
    """
    An identity workflow to check how ideal inverse transform
    affects final evaluation scores.
    """
    wf = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
                        fields=['in_orig', 'in_dist', 'in_tpms', 'in_surf',
                                'in_mask', 'in_field', 'grid_size']),
                        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
                         fields=['out_corr', 'out_tpms',
                                 'out_surf', 'out_field', 'out_mask']),
                         name='outputnode')

    # Invert field
    inv = pe.Node(InverseField(), name='InvertField')

    # Compute corrected images
    merge = pe.Node(niu.Merge(2), name='Merge')
    split = pe.Node(niu.Split(splits=[2, n_tissues]), name='Split')

    # Apply tfm to tpms
    applytfm = pe.Node(FieldBasedWarp(), name="ApplyWarp")

    # Connect
    wf.connect([
        (inputnode,       inv, [('in_field', 'in_field')]),
        (inputnode,     merge, [('in_dist', 'in1'),
                                ('in_tpms', 'in2')]),
        (inputnode,  applytfm, [('in_mask', 'in_mask'),
                                ('in_surf', 'in_surf'),
                                ('grid_size', 'grid_size')]),
        (merge,      applytfm, [('out', 'in_file')]),
        # (inputnode,  applytfm, [('in_field', 'in_field')]),
        (inv,        applytfm, [('out_field', 'in_field')]),
        (applytfm,      split, [('out_file', 'inlist')]),
        (split,    outputnode, [('out1', 'out_corr'),
                                ('out2', 'out_tpms')]),
        (inv,      outputnode, [('out_field', 'out_field')]),
        (applytfm, outputnode, [('out_surf', 'out_surf'),
                                ('out_mask', 'out_mask')])
    ])

    return wf
