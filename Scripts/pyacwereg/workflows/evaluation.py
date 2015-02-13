#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: Oscar Esteban - code@oscaresteban.es
# @Date:   2014-03-12 16:59:14
# @Last Modified by:   oesteban
# @Last Modified time: 2015-02-13 15:02:55

import os
import os.path as op
import numpy as np

import nipype.pipeline.engine as pe             # pipeline engine
from nipype.interfaces import io as nio              # Data i/o
from nipype.interfaces import utility as niu         # utility
from nipype.algorithms import misc as namisc         # misc algorithms
from nipype.algorithms.misc import NormalizeProbabilityMapSet as Normalize
from nipype.algorithms import mesh as namesh
from nipype.algorithms import metrics as namev
from nipype.interfaces import fsl as fsl
from nipype.interfaces import freesurfer as fs

from pyacwereg.interfaces.warps import InverseField
from pyacwereg.interfaces.utility import (ExportSlices, HausdorffDistance,
                                          ComputeEnergy)
from pyacwereg.workflows.model import generate_phantom
from registration import identity_wf, default_regseg


def bspline(name='BSplineEvaluation', shapes=['gyrus'], snr_list=[300],
            N=1, methods=None, results=None):
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
        outputnode.out_surf - the original priors after distortion
          (if available)
        outputnode.out_disp - the displacement field, at image grid resoluton

    """
    wf = pe.Workflow(name=name)

    if methods is None:
        # methods = [identity_wf(n_tissues=2), default_regseg()]
        methods = [default_regseg()]
    else:
        methods = np.atleast_1d(methods).tolist()

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['grid_size', 'out_csv', 'lo_matrix', 'rep_id',
                'hi_matrix', 'snr', 'cortex', 'shape']),
        name='inputnode')

    shapes = np.atleast_1d(shapes).tolist()
    inputnode.iterables = [('shape', shapes),
                           ('snr', snr_list),
                           ('rep_id', range(N))]

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_tpms', 'out_surfs', 'out_field', 'out_coeff',
                'out_overlap']), name='outputnode')

    phantom = generate_phantom()
    wf.connect([
        (inputnode,  phantom, [('shape', 'inputnode.shape'),
                               ('grid_size', 'inputnode.grid_size'),
                               ('lo_matrix', 'inputnode.lo_matrix'),
                               ('hi_matrix', 'inputnode.hi_matrix'),
                               ('snr', 'inputnode.snr'),
                               ('cortex', 'inputnode.cortex'),
                               ('rep_id', 'inputnode.repetition_id')])
    ])

    regseg_low = default_regseg('REGSEG_low')
    ev_regseg_low = registration_ev(name=('Ev_%s' % regseg_low.name))
    ev_regseg_low.inputs.infonode.method = 'REGSEG'
    ev_regseg_low.inputs.infonode.resolution = 'lo'
    norm_low = pe.Node(Normalize(), name='NormalizeFinal_low')
    export0 = pe.Node(ExportSlices(all_axis=True), name='Export_lo')
    sel0 = pe.Node(niu.Select(index=[0]), name='SelectT1w_lo')

    wf.connect([
        (inputnode, ev_regseg_low, [
            ('shape', 'infonode.shape'),
            ('snr', 'infonode.snr'),
            ('rep_id', 'infonode.repetition')]),
        (phantom, ev_regseg_low, [
            ('refnode.out_signal',    'refnode.in_imag'),
            ('refnode.out_tpms',    'refnode.in_tpms'),
            ('out_lowres.out_surfs',   'refnode.in_surf'),
            ('refnode.out_mask',    'refnode.in_mask'),
            ('out_lowres.out_field', 'refnode.in_field')]),
        (phantom, regseg_low, [
            ('refnode.out_surfs', 'inputnode.in_surf'),
            ('out_lowres.out_signal', 'inputnode.in_fixed'),
            ('out_lowres.out_tpms', 'inputnode.in_tpms'),
            ('out_lowres.out_mask', 'inputnode.in_mask')]),
        (regseg_low, norm_low, [
            ('outputnode.out_tpms', 'in_files')]),
        (regseg_low, ev_regseg_low, [
            ('outputnode.out_corr', 'tstnode.in_imag'),
            ('outputnode.out_surf', 'tstnode.in_surf'),
            ('outputnode.out_field', 'tstnode.in_field')]),
        (norm_low, ev_regseg_low, [
            ('out_files', 'tstnode.in_tpms')]),
        (phantom, sel0, [
            ('out_lowres.out_signal', 'inlist')]),
        (sel0, export0, [
            ('out', 'reference')]),
        (phantom, export0, [
            ('out_lowres.out_surfs', 'surfaces0')]),
        (regseg_low, export0, [
            ('outputnode.out_surf', 'surfaces1')])
    ])

    # Connect results output file
    if results is not None:
        ev_regseg_low.inputs.infonode.out_csv = results
    else:
        wf.connect([
            (inputnode, ev_regseg_low, [('out_csv', 'infonode.out_csv')])
        ])

    regseg_hi = default_regseg('REGSEG_hi')
    ev_regseg_hi = registration_ev(name=('Ev_%s' % regseg_hi.name))
    ev_regseg_hi.inputs.infonode.method = 'REGSEG'
    ev_regseg_hi.inputs.infonode.resolution = 'hi'
    norm_hi = pe.Node(Normalize(), name='NormalizeFinal_hi')
    export1 = pe.Node(ExportSlices(all_axis=True), name='Export_hi')
    sel1 = pe.Node(niu.Select(index=[0]), name='SelectT1w_hi')

    wf.connect([
        (inputnode, ev_regseg_hi, [
            ('shape', 'infonode.shape'),
            ('snr', 'infonode.snr'),
            ('rep_id', 'infonode.repetition')]),
        (phantom, ev_regseg_hi, [
            ('refnode.out_signal',    'refnode.in_imag'),
            ('refnode.out_tpms',    'refnode.in_tpms'),
            ('out_hires.out_surfs',   'refnode.in_surf'),
            ('refnode.out_mask',    'refnode.in_mask'),
            ('out_hires.out_field', 'refnode.in_field')]),
        (phantom, regseg_hi, [
            ('refnode.out_surfs', 'inputnode.in_surf'),
            ('out_hires.out_signal', 'inputnode.in_fixed'),
            ('out_hires.out_tpms', 'inputnode.in_tpms'),
            ('out_hires.out_mask', 'inputnode.in_mask')]),
        (regseg_hi, norm_hi, [
            ('outputnode.out_tpms', 'in_files')]),
        (regseg_hi, ev_regseg_hi, [
            ('outputnode.out_corr', 'tstnode.in_imag'),
            ('outputnode.out_surf', 'tstnode.in_surf'),
            ('outputnode.out_field', 'tstnode.in_field')]),
        (norm_hi, ev_regseg_hi, [
            ('out_files', 'tstnode.in_tpms')]),
        (phantom, sel1, [
            ('out_hires.out_signal', 'inlist')]),
        (sel0, export1, [
            ('out', 'reference')]),
        (phantom, export1, [
            ('out_hires.out_surfs', 'surfaces0')]),
        (regseg_hi, export1, [
            ('outputnode.out_surf', 'surfaces1')])
    ])

    # Connect results output file
    if results is not None:
        ev_regseg_hi.inputs.infonode.out_csv = results
    else:
        wf.connect([
            (inputnode, ev_regseg_hi, [('out_csv', 'infonode.out_csv')])
        ])

    # Connect in_field in case it is an identity workflow
    # if 'in_field' in [item[0] for item in reg.inputs.inputnode.items()]:
    #     wf.connect(phantom, 'out_lowres.out_field',
    #                reg, 'inputnode.in_field')

    return wf


def registration_ev(name='EvaluateMapping'):
    """
    Workflow that provides different scores comparing two registration methods.
    It compares images similarity, displacement fields difference,
    mesh distances, and overlap indices.
    """

    def _stats(in_file):
        import numpy as np
        import nibabel as nb
        data = nb.load(in_file).get_data()
        if (np.all(data < 1.0e-5)):
            return [0.0] * 5
        data = np.ma.masked_equal(data, 0)
        result = np.array([data.mean(), data.std(), data.max(), data.min(),
                           np.ma.extras.median(data)])
        return result.tolist()

    def _get_id(inlist):
        return range(len(inlist))

    input_ref = pe.Node(niu.IdentityInterface(
        fields=['in_imag', 'in_tpms', 'in_surf', 'in_field', 'in_mask']),
        name='refnode')
    input_tst = pe.Node(niu.IdentityInterface(
        fields=['in_imag', 'in_tpms', 'in_surf', 'in_field']),
        name='tstnode')
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['snr', 'shape', 'method', 'repetition',
                'resolution', 'out_csv']),
        name='infonode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_tpm_diff', 'out_field_err']),
        name='outputnode')
    merge_ref = pe.Node(fsl.Merge(dimension='t'), name='ConcatRefInputs')
    merge_tst = pe.Node(fsl.Merge(dimension='t'), name='ConcatTestInputs')
    overlap = pe.Node(namev.FuzzyOverlap(weighting='volume'), name='Overlap')
    diff_im = pe.Node(namev.Similarity(metric='cc'), name='ContrastDiff')
    inv_fld = pe.Node(InverseField(), name='InvertField')
    diff_fld = pe.Node(namev.ErrorMap(), name='FieldDiff')
    mesh = pe.MapNode(HausdorffDistance(cells_mode=True),
                      iterfield=['surface1', 'surface2'],
                      name='SurfDistance')
    csv = pe.MapNode(namisc.AddCSVRow(infields=['surf_id', 'surfdist_avg']),
                     name="AddRow", iterfield=['surf_id', 'surfdist_avg'])
    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,        csv, [('shape', 'model_type'),
                                 ('snr', 'snr'),
                                 ('method', 'method'),
                                 ('resolution', 'resolution'),
                                 ('repetition', 'repetition'),
                                 ('out_csv', 'in_file')]),
        # (input_ref,  merge_ref, [('in_imag', 'in_files')]),
        # (input_tst,  merge_tst, [('in_imag', 'in_files')]),
        # (input_ref,    overlap, [('in_tpms', 'in_ref')]),
        # (input_tst,    overlap, [('in_tpms', 'in_tst')]),
        # (input_ref,    diff_im, [('in_mask', 'mask1'),
        #                          ('in_mask', 'mask2')]),
        # (merge_ref,    diff_im, [('merged_file', 'volume1')]),
        # (merge_tst,    diff_im, [('merged_file', 'volume2')]),
        # (input_ref,    inv_fld, [('in_field', 'in_field')]),
        # (input_ref,   diff_fld, [('in_mask', 'mask')]),
        # (inv_fld,     diff_fld, [('out_field', 'in_ref')]),
        # (input_tst,   diff_fld, [('in_field', 'in_tst')]),
        # (overlap,          csv, [('jaccard', 'fji_avg'),
        #                          ('class_fji', 'fji_tpm'),
        #                          ('dice', 'fdi_avg'),
        #                          ('class_fdi', 'fdi_tpm')]),
        # (diff_im,          csv, [('similarity', 'cc_image')]),
        # (diff_fld,         csv, [(('out_map', _stats), 'fmap_error')]),
        # (csv,       outputnode, [('csv_file', 'out_file')]),
        # (overlap,   outputnode, [('diff_file', 'out_tpm_diff')]),
        # (diff_fld,  outputnode, [('out_map', 'out_field_err')]),
        (input_ref,       mesh, [('in_surf', 'surface1')]),
        (input_tst,       mesh, [('in_surf', 'surface2')]),
        (mesh,             csv, [('avg_hd', 'surfdist_avg'),
                                 (('avg_hd', _get_id), 'surf_id')])
        # (mesh,             csv, [('max_hd', 'surfdist_hausdorff'),
        #                          ('avg_hd', 'surfdist_avg'),
        #                          ('std_hd', 'surfdist_std'),
        #                          ('stats_hd', 'surfdist_stats')])
    ])
    return wf


def map_energy(name='EnergyMapping'):

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['reference', 'surfaces0', 'surfaces1', 'in_mask']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['desc_zero', 'out_diff']), name='outputnode')

    ref_e = pe.Node(ComputeEnergy(), name='ComputeZeroEnergy')
    diff = pe.MapNode(namesh.ComputeMeshWarp(), name='ComputeError',
                      iterfield=['surface1', 'surface2'])

    mapper = warp_n_map()
    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,     ref_e,  [('reference', 'reference'),
                                 ('surfaces0', 'surfaces'),
                                 ('in_mask', 'in_mask')]),
        (ref_e,     outputnode, [('out_file', 'desc_zero')]),
        (inputnode,       diff, [('surfaces0', 'surface1'),
                                 ('surfaces1', 'surface2')]),
        (diff,      outputnode, [('out_warp', 'out_diff')]),

        (inputnode,     mapper, [('reference', 'inputnode.reference'),
                                 ('in_mask', 'inputnode.in_mask')]),
        (diff,          mapper, [('out_warp', 'inputnode.surf_warp')]),
        (ref_e,         mapper, [('out_desc', 'inputnode.descriptors')])
    ])
    return wf


def warp_n_map(name='EnergyWarpAndMap', out_csv='energies.csv'):
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['reference', 'surf_warp', 'in_mask', 'errfactor',
                'descriptors']), name='inputnode')
    inputnode.iterables = ('errfactor', [-0.05, 0.05])

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_energy']), name='outputnode')

    applyef = pe.MapNode(
        namesh.MeshWarpMaths(operation='mul'), name='MeshMaths',
        iterfield=['in_surf'])
    mapeneg = pe.Node(ComputeEnergy(), name='ComputeEnergy')
    getval = pe.Node(nio.JSONFileGrabber(), name='GetEnergy')

    csv = pe.Node(namisc.AddCSVRow(in_file=out_csv),
        name="AddRow")

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,    applyef, [('surf_warp', 'in_surf'),
                                 ('errfactor', 'operator')]),
        (applyef,      mapeneg, [('out_file', 'surfaces')]),
        (inputnode,    mapeneg, [('reference', 'reference'),
                                 ('in_mask', 'in_mask'),
                                 ('descriptors', 'descriptors')]),
        (mapeneg,       getval, [('out_file', 'in_file')]),
        (getval,           csv, [('total', 'total')]),
        (inputnode,        csv, [('errfactor', 'error')]),
        (mapeneg,   outputnode, [('out_file', 'out_energy')])
    ])
    return wf
