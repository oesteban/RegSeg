#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: Oscar Esteban - code@oscaresteban.es
# @Date:   2014-03-12 16:59:14
# @Last Modified by:   oesteban
# @Last Modified time: 2014-04-11 15:19:09

import os
import os.path as op
import numpy as np

import nipype.interfaces.io as nio              # Data i/o
import nipype.interfaces.utility as niu         # utility
import nipype.algorithms.misc as namisc         # misc algorithms
import nipype.algorithms.mesh as namesh
import nipype.algorithms.eval as namev
import nipype.interfaces.fsl as fsl
from nipype.interfaces.nipy.utils import Similarity
import nipype.pipeline.engine as pe             # pipeline engine

import pyacwereg.nipype.interfaces as iface
from pyacwereg.utils.misc import normalize_tpms as normalize

from smri import prepare_smri
from distortion import bspline_deform
from registration import identity_wf,default_regseg

def registration_ev( name='EvaluateMapping' ):
    """ Workflow that provides different scores comparing two registration methods. It compares images
    similarity, displacement fields difference, mesh distances, and overlap indices.
    """

    def _stats( in_file ):
        import numpy as np
        import nibabel as nb

        data = nb.load( in_file ).get_data()

        if ( np.all( data<1.0e-5 ) ):
          return [ 0.0 ]*5
        data = np.ma.masked_equal( data, 0 )
        result = np.array([ data.mean(), data.std(), data.max(), data.min(), np.ma.extras.median(data)])
        return result.tolist()


    wf = pe.Workflow( name=name )
    input_ref = pe.Node( niu.IdentityInterface( fields=[ 'in_imag',
                        'in_tpms','in_surf','in_field', 'in_mask' ] ),
                        name='refnode' )
    input_tst = pe.Node( niu.IdentityInterface( fields=['in_imag', 'in_tpms',
                                                        'in_surf','in_field' ] ),
                        name='tstnode' )
    inputnode = pe.Node( niu.IdentityInterface( fields=['subject_id', 'method', 'out_csv']),
                         name='infonode' )
    outputnode = pe.Node(niu.IdentityInterface(fields=[ 'out_file', 'out_tpm_diff', 'out_field_err' ]),
                         name='outputnode' )

    merge_ref = pe.Node( fsl.Merge(dimension='t'), name='ConcatRefInputs' )
    merge_tst = pe.Node( fsl.Merge(dimension='t'), name='ConcatTestInputs' )


    overlap = pe.Node( namev.FuzzyOverlap(weighting='volume'), name='Overlap' )
    diff_im = pe.Node( Similarity(metric='cc'), name='ContrastDiff')
    inv_fld = pe.Node( iface.InverseField(), name='InvertField' )
    diff_fld = pe.Node( namev.ErrorMap(), name='FieldDiff')
#    mesh = pe.MapNode( namesh.P2PDistance(weighting='surface'),
#                      iterfield=[ 'surface1','surface2' ],
#                      name='SurfDistance')

    csv = pe.Node( namisc.AddCSVRow(), name="AddRow" )

    wf.connect( [
                ( inputnode,         csv, [( 'subject_id', 'subject_id'), ('method','method'),
                                           ( 'out_csv', 'in_file' )])
               ,( input_ref,   merge_ref, [( 'in_imag', 'in_files' )])
               ,( input_tst,   merge_tst, [( 'in_imag', 'in_files' )])
               ,( input_ref,     overlap, [( 'in_tpms', 'in_ref')] )
               ,( input_tst,     overlap, [( 'in_tpms', 'in_tst')] )
               ,( input_ref,     diff_im, [( 'in_mask','mask1'),
                                           ( 'in_mask','mask2')])
               ,( merge_ref,     diff_im, [( 'merged_file', 'volume1')])
               ,( merge_tst,     diff_im, [( 'merged_file', 'volume2')])
               ,( input_ref,     inv_fld, [( 'in_field', 'in_field') ])
               ,( input_ref,    diff_fld, [( 'in_mask','mask')])
               ,( inv_fld,      diff_fld, [( 'out_field', 'in_ref') ])
               ,( input_tst,    diff_fld, [( 'in_field', 'in_tst')])
#               ,( input_ref,        mesh, [( 'in_surf', 'surface1')])
#               ,( input_tst,        mesh, [( 'in_surf', 'surface2')])
               ,( overlap,           csv, [( 'jaccard', 'fji_avg'), ('class_fji','fji_tpm'),
                                           ( 'dice', 'fdi_avg'), ('class_fdi', 'fdi_tpm') ])
               ,( diff_im,           csv, [( 'similarity','cc_image')])
               ,( diff_fld,          csv, [(('out_map',_stats), 'fmap_error' ) ])
#               ,( mesh,              csv, [( 'distance', 'surf_dist')])
               ,( csv,        outputnode, [( 'csv_file', 'out_file')])
               ,( overlap,    outputnode, [( 'diff_file','out_tpm_diff')])
               ,( diff_fld,   outputnode, [( 'out_map', 'out_field_err')])
    ])

    return wf

def bspline( name='BSplineEvaluation', methods=None, results=None ):
    """ A workflow to evaluate registration methods generating a gold standard
    with random bspline deformations.

    A list of nipype workflows can be plugged-in, using the methods input. If
    methods is None, then a default regseg method is run.


    Inputs in methods workflows
    ---------------------------

    methods workflows must define the following inputs:
        inputnode.in_surf - the input prior / surfaces in orig space
        inputnode.in_dist - the distorted images
        inputnode.in_tpms - the distorted TPMs (tissue probability maps)
        inputnode.in_orig - the original images, undistorted


    Outputs in methods workflows
    ----------------------------

        outputnode.out_corr - the distorted images, after correction
        outputnode.out_tpms - the corrected TPMs
        outputnode.out_surf - the original priors after distortion (if available)
        outputnode.out_disp - the displacement field, at image grid resoluton

    """
    wf = pe.Workflow( name=name )


    if methods is None:
        methods = [ identity_wf(), default_regseg() ]
    else:
        methods = np.atleast_1d( methods ).tolist()

    inputnode = pe.Node( niu.IdentityInterface( fields=[ 'subject_id', 'data_dir',
                                                         'grid_size', 'out_csv' ] ),
                        name='inputnode' )
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file', 'out_tpms',
                         'out_surfs','out_field', 'out_coeff', 'out_overlap' ]),
                         name='outputnode' )

    prep = prepare_smri()
    dist = bspline_deform()

    wf.connect([
             ( inputnode,  prep, [ ('subject_id','inputnode.subject_id'),('data_dir','inputnode.data_dir') ])
            ,( inputnode,  dist, [ ('grid_size', 'inputnode.grid_size')])
            ,( prep,       dist, [ ('outputnode.out_smri_brain','inputnode.in_file'),
                                   ('outputnode.out_surfs', 'inputnode.in_surfs'),
                                   ('outputnode.out_tpms', 'inputnode.in_tpms')])
            ,( dist, outputnode, [ ('outputnode.out_file','out_file'),
                                   ('outputnode.out_field','out_field'),
                                   ('outputnode.out_coeff','out_coeff')])
    ])

    evwfs = []
    norm_tpms = []
    for i,reg in enumerate(methods):
        evwfs.append( registration_ev( name=('Ev_%s' % reg.name) ) )
        evwfs[i].inputs.infonode.method = reg.name

        norm_tpms.append(pe.Node( niu.Function( input_names=['in_files','in_mask'],
                         output_names=['out_files'], function=normalize ),
                         name='Normalize%02d' % i ) )

        wf.connect( [
             (inputnode,    evwfs[i], [('subject_id', 'infonode.subject_id') ] )
            ,(prep,              reg, [('outputnode.out_surfs','inputnode.in_surf'),
                                       ('outputnode.out_smri_brain', 'inputnode.in_orig' ),
                                       ('outputnode.out_mask', 'inputnode.in_mask' ) ])
            ,(dist,              reg, [('outputnode.out_file', 'inputnode.in_dist'),
                                       ('outputnode.out_tpms', 'inputnode.in_tpms') ])
            ,(prep,         evwfs[i], [('outputnode.out_smri_brain', 'refnode.in_imag'),
                                       ('outputnode.out_tpms',       'refnode.in_tpms'),
                                       ('outputnode.out_surfs',      'refnode.in_surf'),
                                       ('outputnode.out_mask',       'refnode.in_mask'), ])
            ,(dist,         evwfs[i], [('outputnode.out_field',      'refnode.in_field' ) ])
            ,(reg,      norm_tpms[i], [('outputnode.out_tpms', 'in_files') ] )
            ,(reg,          evwfs[i], [('outputnode.out_corr', 'tstnode.in_imag'),
                                       ('outputnode.out_surf', 'tstnode.in_surf'),
                                       ('outputnode.out_field','tstnode.in_field' ) ])
            ,(norm_tpms[i], evwfs[i], [('out_files', 'tstnode.in_tpms')])
        ])

        # Connect in_field in case it is an identity workflow
        if 'in_field' in [ item[0] for item in reg.inputs.inputnode.items() ]:
            wf.connect( [
                ( dist, reg, [( 'outputnode.out_field', 'inputnode.in_field')])
            ])

        # Connect results output file
        if not results is None:
            evwfs[i].inputs.infonode.out_csv = results
        else:
            wf.connect([
             ( inputnode, evwfs[i], [ ('out_csv', 'infonode.out_csv') ])
            ])


    return wf
