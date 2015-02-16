#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban - code@oscaresteban.es
# @Date:   2014-03-28 20:38:30
# @Last Modified by:   Oscar Esteban
# @Last Modified time: 2015-02-16 14:47:13

import os
import os.path as op

import nipype.pipeline.engine as pe             # pipeline engine
from nipype.interfaces import io as nio              # Data i/o
from nipype.interfaces import utility as niu         # utility
from nipype.interfaces import fsl

from pyacwereg.interfaces.acwereg import ACWEReg, ACWEReport
from pyacwereg.interfaces.warps import FieldBasedWarp, InverseField


def regseg_wf(name='REGSEG', enhance_inputs=True, usemask=False):
    wf = pe.Workflow(name=name)

    regseg_inputs = ['iterations', 'alpha', 'beta', 'step_size',
                     'grid_spacing', 'convergence_energy',
                     'convergence_window', 'f_smooth',
                     'images_verbosity', 'scales', 'descript_update',
                     'convergence_value', 'descript_adaptative']

    wf_inputs = ['in_fixed', 'in_tpms', 'in_surf', 'in_mask']
    inputnode = pe.Node(niu.IdentityInterface(
                        fields=regseg_inputs + wf_inputs), name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_corr', 'out_tpms', 'out_enh', 'reg_msk', 'out_surf',
                'out_field', 'out_mask']),
        name='outputnode')

    # Registration
    regseg = pe.Node(ACWEReg(), name="ACWERegistration")
    report = pe.Node(ACWEReport(), name="ACWEReport")

    # Apply tfm to tpms
    applytfm = pe.MapNode(FieldBasedWarp(), name="ApplyWarp",
                          iterfield=['in_file'])

    # Connect
    wf.connect([
        (inputnode,   regseg, [(f, f) for f in regseg_inputs]),
        (inputnode,   regseg, [('in_surf', 'in_prior')]),
        (inputnode, applytfm, [('in_tpms', 'in_file'),
                               ('in_mask', 'in_mask')]),
        (regseg,    applytfm, [('out_field', 'in_field')]),
        (regseg,      report, [('out_log', 'in_log')]),
        (regseg,  outputnode, [('out_warped', 'out_corr'),
                               ('out_field', 'out_field'),
                               ('out_surfs', 'out_surf')]),
        (applytfm, outputnode, [('out_file', 'out_tpms'),
                                ('out_mask', 'out_mask')])
    ])

    if usemask:
        dilate = pe.Node(fsl.maths.MathsCommand(
            nan2zeros=True, args='-kernel sphere 5 -dilM'), name='MskDilate')
        wf.connect([(inputnode,   dilate, [('in_mask', 'in_file')]),
                    (dilate,      regseg, [('out_file', 'in_mask')]),
                    (dilate,  outputnode, [('out_file', 'reg_msk')])])
    else:
        zeromsk = pe.Node(niu.Function(
            input_names=['in_file'], output_names=['out_file'],
            function=_gen_zmsk), name='ZeroMsk')
        wf.connect([(inputnode,   zeromsk, [('in_mask', 'in_file')]),
                    (zeromsk,  outputnode, [('out_file', 'reg_msk')])])

    if enhance_inputs:
        enh = pe.MapNode(niu.Function(
            function=enh_image, input_names=['in_file'],
            output_names=['out_file']), iterfield=['in_file'], name='Enhance')
        wf.connect([
            (inputnode,        enh, [('in_fixed', 'in_file')]),
            (enh,           regseg, [('out_file', 'in_fixed')]),
            (enh,       outputnode, [('out_file', 'out_enh')])
        ])
    else:
        wf.connect([
            (inputnode,     regseg, [('in_fixed', 'in_fixed')]),
            (inputnode, outputnode, [('in_fixed', 'out_enh')])
        ])

    return wf


def default_regseg(name='REGSEGDefault'):
    wf = regseg_wf(name=name, enhance_inputs=False)

    # Registration
    wf.inputs.inputnode.iterations = [500, 250]
    wf.inputs.inputnode.descript_update = [None, None]
    wf.inputs.inputnode.step_size = [1.e-3, 0.01]
    wf.inputs.inputnode.alpha = [0.0, 0.0]
    wf.inputs.inputnode.beta = [0.0, 0.0]
    wf.inputs.inputnode.grid_spacing = [16., 8.]
    wf.inputs.inputnode.convergence_energy = [True] * 2
    wf.inputs.inputnode.descript_adaptative = [True, False]
    wf.inputs.inputnode.convergence_window = [60, 5]
    wf.inputs.inputnode.f_smooth = [None, None]
    wf.inputs.inputnode.images_verbosity = 1
    wf.inputs.inputnode.convergence_value = [1.0e-6, 1.0e-8]
    return wf


def t2b_workflow(name='T2B', minimal=False, enc_dir=None,
                 enhance_b0=True, icorr=True):
    """
    The T2w-registration based method (T2B) implements an SDC by nonlinear
    registration of the anatomically correct *T2w* image to the *b0* image
    of the *dMRI* dataset. The implementation here tries to reproduce the one
    included in ExploreDTI `(Leemans et al., 2009)
    <http://www.exploredti.com/ref/ExploreDTI_ISMRM_2009.pdf>`_, which is
    also based on `(Irfanoglu et al., 2012)
    <http://dx.doi.org/10.1016/j.neuroimage.2012.02.054>`_.

    :param str name: a unique name for the workflow.

    :inputs:

        * in_t2w: the reference T2w image

    :outputs:

        * outputnode.corrected_image: the dMRI image after correction


    Example::

    >>> t2b = t2b_workflow()
    >>> t2b.inputs.inputnode.in_dwi = 'dwi_brain.nii'
    >>> t2b.inputs.inputnode.in_bval = 'dwi.bval'
    >>> t2b.inputs.inputnode.in_mask = 'b0_mask.nii'
    >>> t2b.inputs.inputnode.in_t2w = 't2w_brain.nii'
    >>> t2b.run() # doctest: +SKIP

    """
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_dwi', 'in_bval', 'in_t2w', 'dwi_mask', 't2w_mask',
                'in_param']), name='inputnode')

    avg_b0 = pe.Node(pmisc.ComputeAveragedB0(), name='avg_b0')
    cache_b0 = pe.Node(niu.IdentityInterface(fields=['b0', 'mask']),
                       name='B0Cache')
    enh_t2 = pe.Node(pmisc.SigmoidFilter(upper_perc=92.0, lower_perc=65.0),
                     name='enh_T2')
    reg_param = pe.Node(nii.JSONFileGrabber(defaults={'enc_dir': 'y'}),
                        name='T2BSettings')
    reg = pe.Node(elastix.Registration(num_threads=1), name='Elastix')
    tfx_b0 = pe.Node(elastix.EditTransform(), name='tfm_b0')
    split_dwi = pe.Node(fsl.utils.Split(dimension='t'), name='split_dwi')
    warp = pe.MapNode(elastix.ApplyWarp(), iterfield=['moving_image'],
                      name='UnwarpDWIs')
    warp_prop = pe.Node(elastix.AnalyzeWarp(), name='DisplFieldAnalysis')
    warpbuff = pe.Node(niu.IdentityInterface(fields=['unwarped']),
                       name='UnwarpedCache')
    mskdwis = pe.MapNode(fs.ApplyMask(), iterfield='in_file', name='MaskDWIs')
    thres = pe.MapNode(Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegs')
    merge_dwi = pe.Node(fsl.utils.Merge(dimension='t'), name='merge_dwis')
    tfx_msk = pe.Node(elastix.EditTransform(
        interpolation='nearest', output_type='unsigned char'),
        name='MSKInterpolator')
    corr_msk = pe.Node(elastix.ApplyWarp(), name='UnwarpMsk')
    closmsk = pe.Node(fsl.maths.MathsCommand(
        nan2zeros=True, args='-kernel sphere 3 -dilM -kernel sphere 2 -ero'),
        name='MaskClosing')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,    avg_b0, [('in_dwi', 'in_dwi'),
                                ('in_bval', 'in_bval')]),
        (inputnode,   selbmsk, [('in_segs', 'inlist')]),
        (inputnode, split_dwi, [('in_dwi', 'in_file')]),
        (inputnode,  corr_msk, [('in_mask', 'moving_image')]),
        (inputnode,    enh_t2, [('in_t2w', 'in_file')]),
        (selbmsk,      enh_t2, [('out', 'in_mask')]),
        (reg_param,       reg, [('out_file', 'parameters')]),
        (enh_t2,          reg, [('out_file', 'fixed_image')]),
        (cache_b0,        reg, [('b0', 'moving_image'),
                                ('mask', 'moving_mask')]),
        (selbmsk,         reg, [('out', 'fixed_mask')]),
        (reg,          tfx_b0, [(('transform', _get_last), 'transform_file')]),
        (avg_b0,       tfx_b0, [('out_file', 'reference_image')]),
        (tfx_b0,    warp_prop, [('output_file', 'transform_file')]),
        (tfx_b0,         warp, [('output_file', 'transform_file')]),
        (split_dwi,      warp, [('out_files', 'moving_image')]),
        (warpbuff,    mskdwis, [('unwarped', 'in_file')]),
        (closmsk,     mskdwis, [('out_file', 'mask_file')]),
        (mskdwis,       thres, [('out_file', 'in_file')]),
        (thres,     merge_dwi, [('out_file', 'in_files')]),
        (reg,         tfx_msk, [(('transform', _get_last), 'transform_file')]),
        (avg_b0,      tfx_msk, [('out_file', 'reference_image')]),
        (tfx_msk,    corr_msk, [('output_file', 'transform_file')]),
        (corr_msk,    closmsk, [('warped_file', 'in_file')]),
        (merge_dwi,     'out', [('merged_file', 'dwi')]),
        (closmsk,       'out', [('out_file', 'dwi_mask')]),
        (warp_prop,     'out', [('jacdet_map', 'jacobian')])
    ])

    if icorr:
        jac_mask = pe.Node(fs.ApplyMask(), name='mask_jac')
        mult = pe.MapNode(MultiImageMaths(op_string='-mul %s'),
                          iterfield=['in_file'], name='ModulateDWIs')
        wf.connect([
            (closmsk,      jac_mask, [('out_file', 'mask_file')]),
            (warp_prop,    jac_mask, [('jacdet_map', 'in_file')]),
            (warp,             mult, [('warped_file', 'in_file')]),
            (jac_mask,         mult, [('out_file', 'operand_files')]),
            (mult,         warpbuff, [('out_file', 'unwarped')])
        ])
    else:
        wf.connect([
            (warp,         warpbuff, [('warped_file', 'unwarped')])
        ])

    if enhance_b0:
        enh_b0 = pe.Node(pmisc.EnhanceB0(), name='enh_b0')
        wf.connect([
            (inputnode,      enh_b0, [('in_mask', 'in_mask')]),
            (avg_b0,         enh_b0, [('out_file', 'in_file')]),
            (enh_b0,       cache_b0, [('out_file', 'b0'),
                                      ('out_mask', 'mask')])
        ])
    else:
        wf.connect([
            (inputnode,    cache_b0, [('in_mask', 'mask')]),
            (avg_b0,       cache_b0, [('out_file', 'b0')])
        ])

    if not minimal:
        invwarp = pe.Node(fsl.InvWarp(relative=True), name='InvWarp')
        warpsrf = pe.MapNode(WarpPoints(), iterfield=['points'],
                             name='UnwarpSurfs')

        tfx_tpms = pe.Node(elastix.EditTransform(interpolation='linear',
                                                 output_type='float'),
                           name='TPMInterpolator')
        corr_tpm = pe.MapNode(elastix.ApplyWarp(), iterfield=['moving_image'],
                              name='UnwarpTPMs')
        normtpm = pe.Node(NormalizeTPM(), name='NormTPMs')
        tfx_segs = pe.Node(elastix.EditTransform(interpolation='nearest',
                                                 output_type='unsigned short'),
                           name='NNInterpolator')
        corr_segs = pe.MapNode(elastix.ApplyWarp(), iterfield=['moving_image'],
                               name='UnwarpSegs')

        wf.connect([
            (reg,       tfx_segs, [(
                ('transform', _get_last), 'transform_file')]),
            (tfx_segs, corr_segs, [('output_file', 'transform_file')]),
            (inputnode, corr_segs, [('in_segs', 'moving_image')]),
            (inputnode,  warpsrf, [('in_surf', 'points')]),
            (warp_prop,  invwarp, [('disp_field', 'warp')]),
            (avg_b0,     invwarp, [('out_file', 'reference')]),
            (invwarp,    warpsrf, [('inverse_warp', 'warp')]),
            (reg,       tfx_tpms, [
                (('transform', _get_last), 'transform_file')]),
            (tfx_tpms,  corr_tpm, [('output_file', 'transform_file')]),
            (inputnode, corr_tpm, [('in_tpms', 'moving_image')]),
            (corr_segs,    'out', [('warped_file', 'segs')]),
            (corr_tpm,   normtpm, [('warped_file', 'in_files')]),
            (selbmsk,    normtpm, [('out', 'in_mask')]),
            (normtpm,      'out', [('out_files', 'tpms')]),
            (warpsrf,      'out', [('out_points', 'surf')])
        ])

    return wf


def identity_wf(name='Identity', n_tissues=3):
    """
    An identity workflow to check how ideal inverse transform
    affects final evaluation scores.
    """
    wf = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
                        fields=['in_fixed', 'in_tpms', 'in_surf',
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
        (inputnode,     merge, [('in_fixed', 'in1'),
                                ('in_tpms', 'in2')]),
        (inputnode,  applytfm, [('in_mask', 'in_mask'),
                                ('in_surf', 'in_surf'),
                                ('grid_size', 'grid_size')]),
        (merge,      applytfm, [('out', 'in_file')]),
        (inv,        applytfm, [('out_field', 'in_field')]),
        (applytfm,      split, [('out_file', 'inlist')]),
        (split,    outputnode, [('out1', 'out_corr'),
                                ('out2', 'out_tpms')]),
        (inv,      outputnode, [('out_field', 'out_field')]),
        (applytfm, outputnode, [('out_surf', 'out_surf'),
                                ('out_mask', 'out_mask')])
    ])

    return wf


def enh_image(in_file, irange=2000., out_file=None):
    import numpy as np
    import nibabel as nb
    import os.path as op

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_enh.nii.gz' % fname)

    nii = nb.load(in_file)
    data = nii.get_data()
    data[data < 0] = 0.0
    data[data > 1.0] = 0.0
    imax = data.max()
    data = (irange / imax) * data

    nb.Nifti1Image(data, nii.get_affine(), nii.get_header()).to_filename(
        out_file)
    return out_file


def _gen_zmsk(in_file, out_file=None):
    import numpy as np
    import nibabel as nb
    import os.path as op

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_zmsk.nii.gz' % fname)

    nii = nb.load(in_file)
    msk = np.ones(nii.get_shape()[:3])
    hdr = nii.get_header().copy()
    hdr.set_data_shape(msk.shape)
    hdr.set_data_dtype(np.uint8)
    hdr.set_xyzt_units('mm')

    nb.Nifti1Image(msk.astype(np.uint8), nii.get_affine(), hdr).to_filename(
        out_file)
    return out_file
