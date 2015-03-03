#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2014-10-23 14:43:23
# @Last Modified by:   oesteban
# @Last Modified time: 2015-03-03 15:02:29

import os
import os.path as op

from nipype.pipeline import engine as pe             # pipeline engine
from nipype.interfaces import io as nio              # Data i/o
from nipype.interfaces import utility as niu         # utility
from nipype.interfaces import fsl
from nipype.interfaces import freesurfer as fs
from nipype.algorithms.misc import NormalizeProbabilityMapSet as Normalize

from pysdcev.workflows.distortion import bspline_deform

from pyacwereg.workflows.surfaces import extract_surface
from pyacwereg.interfaces import phantoms as pip
from pyacwereg.interfaces.warps import FieldBasedWarp, InverseField
from pyacwereg.interfaces.utility import Surf2Vol


def generate_phantom(name='PhantomGeneration'):
    """
    A phantom generation workflow
    """
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['shape', 'hi_matrix', 'lo_matrix', 'snr', 'cortex',
                'grid_size', 'repetition_id']),
        name='inputnode')

    out_lowres = pe.Node(niu.IdentityInterface(
        fields=['out_signal', 'out_mask', 'out_tpms', 'out_surfs',
                'out_field', 'out_coeff', 'grid_size']),
        name='out_lowres')

    out_hires = pe.Node(niu.IdentityInterface(
        fields=['out_signal', 'out_mask', 'out_tpms', 'out_surfs',
                'out_field', 'out_coeff', 'grid_size']),
        name='out_hires')

    refnode = pe.Node(niu.IdentityInterface(
        fields=['out_signal', 'out_mask', 'out_tpms', 'out_surfs']),
        name='refnode')

    model = pe.Node(pip.Phantom(), name='GenerateModel')
    split = pe.Node(fsl.Split(dimension='t'), name='Split')
    selm0 = pe.Node(niu.Split(splits=[1, 2], squeeze=True),
                    name='SepModel0')
    selm1 = pe.Node(niu.Split(splits=[1, 1, 1], squeeze=True),
                    name='SepModel1')
    signal0 = pe.Node(pip.SimulateSMRI(), name='Simulate0')
    merge0 = pe.Node(niu.Merge(2), name='SimMerge0')

    surf0 = extract_surface('GenSurf0')
    surf0.inputs.inputnode.labels = [1]
    surf0.inputs.inputnode.name = '00.white'
    surf1 = extract_surface('GenSurf1')
    surf1.inputs.inputnode.labels = [1]
    surf1.inputs.inputnode.name = '01.pial'
    msurf = pe.Node(niu.Merge(2), name='MergeSurfs')

    down = pe.Node(fs.MRIConvert(), name='Downsample')

    dist = bspline_deform(n_tissues=0)

    surf2vol0 = pe.Node(Surf2Vol(), name='Surf2Volume_HR')
    surf2vol1 = pe.Node(Surf2Vol(), name='Surf2Volume_LR')
    norm0 = pe.Node(Normalize(), name='NormalizeTPMs_HR')
    norm1 = pe.Node(Normalize(), name='NormalizeTPMs_LR')

    tpmmsk0 = pe.Node(niu.Split(splits=[2, 1, 1]), name='TPMsSplit_HR')
    tpmmsk1 = pe.Node(niu.Split(splits=[2, 1, 1]), name='TPMsSplit_LR')

    msk0 = pe.Node(niu.Function(function=_bin_n_msk, input_names=['in_files'],
                                output_names=['out_file']), name='binNmsk_HR')
    msk1 = pe.Node(niu.Function(function=_bin_n_msk, input_names=['in_files'],
                                output_names=['out_file']), name='binNmsk_LR')

    selt0 = pe.Node(niu.Split(splits=[1, 1, 1, 1], squeeze=True),
                    name='SeparateTissue_HR')
    selt1 = pe.Node(niu.Split(splits=[1, 1, 1, 1], squeeze=True),
                    name='SeparateTissue_LR')

    merge1 = pe.Node(niu.Merge(2), name='SimMerge_HR')
    merge2 = pe.Node(niu.Merge(2), name='SimMerge_LR')

    signal1 = pe.Node(pip.SimulateSMRI(), name='SimulateHR')
    signal2 = pe.Node(pip.SimulateSMRI(), name='SimulateLR')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,   model,       [('shape', 'shape'),
                                    ('hi_matrix', 'matrix_size'),
                                    ('cortex', 'cortex'),
                                    ('repetition_id', 'seed')]),
        (model,       split,       [('out_file', 'in_file')]),
        (split,       selm1,       [('out_files', 'inlist')]),
        (selm1,       signal0,     [('out1', 'frac_csf'),
                                    ('out2', 'frac_wm'),
                                    ('out3', 'frac_gm')]),
        (signal0,     surf0,       [('out_t1w', 'inputnode.norm')]),
        (selm1,       surf0,       [('out2', 'inputnode.aseg')]),
        (signal0,     surf1,       [('out_t1w', 'inputnode.norm')]),
        (model,       surf1,       [('out_mask', 'inputnode.aseg')]),
        (surf0,       msurf,       [('outputnode.out_surf', 'in1')]),
        (surf1,       msurf,       [('outputnode.out_surf', 'in2')]),
        (split,       selm0,       [('out_files', 'inlist')]),
        (inputnode,   dist,        [('grid_size', 'inputnode.grid_size')]),
        (msurf,       dist,        [('out', 'inputnode.in_surfs')]),
        (model,       dist,        [('out_mask', 'inputnode.in_mask')]),
        (selm0,       dist,        [('out2', 'inputnode.in_file')]),

        (signal0,     surf2vol0,   [('out_t1w', 'reference')]),
        (dist,        surf2vol0,   [('outputnode.out_surfs', 'surfaces')]),
        (surf2vol0,   norm0,       [('out_tpm', 'in_files')]),
        (norm0,       selt0,       [('out_files', 'inlist')]),
        (selt0,       signal1,     [('out1', 'frac_wm'),
                                    ('out2', 'frac_gm'),
                                    ('out3', 'frac_csf')]),
        (inputnode,   signal1,     [('snr', 'snr')]),
        (signal1,     merge1,      [('out_t1w', 'in1'),
                                    ('out_t2w', 'in2')]),
        (norm0,       tpmmsk0,     [('out_files', 'inlist')]),
        (tpmmsk0,     msk0,        [('out1', 'in_files')]),

        (signal0,     down,        [('out_t1w', 'in_file'),
                                    (('out_t1w', _half_voxsize), 'vox_size')]),
        (down,        surf2vol1,   [('out_file', 'reference')]),
        (dist,        surf2vol1,   [('outputnode.out_surfs', 'surfaces')]),
        (surf2vol1,   norm1,       [('out_tpm', 'in_files')]),
        (norm1,       selt1,       [('out_files', 'inlist')]),
        (selt1,       signal2,     [('out1', 'frac_wm'),
                                    ('out2', 'frac_gm'),
                                    ('out3', 'frac_csf')]),
        (inputnode,   signal2,     [('snr', 'snr')]),
        (signal2,     merge2,      [('out_t1w', 'in1'),
                                    ('out_t2w', 'in2')]),
        (norm1,       tpmmsk1,     [('out_files', 'inlist')]),
        (tpmmsk1,     msk1,        [('out1', 'in_files')]),

        # reference outputs
        (signal0,     merge0,      [('out_t1w', 'in1'),
                                    ('out_t2w', 'in2')]),
        (msurf,       refnode,     [('out', 'out_surfs')]),
        (selt0,       refnode,     [('out2', 'out_tpms')]),
        (model,       refnode,     [('out_mask', 'out_mask')]),
        (merge0,      refnode,     [('out', 'out_signal')]),

        # distorted outputs
        (inputnode,   out_hires,   [('grid_size', 'grid_size')]),
        (merge1,      out_hires,   [('out', 'out_signal')]),
        (msk0,        out_hires,   [('out_file', 'out_mask')]),
        (tpmmsk0,     out_hires,   [('out1', 'out_tpms')]),
        (dist,        out_hires,   [('outputnode.out_field', 'out_field'),
                                    ('outputnode.out_coeff', 'out_coeff'),
                                    ('outputnode.out_surfs', 'out_surfs')]),

        # distorted outputs
        (inputnode,   out_lowres,  [('grid_size', 'grid_size')]),
        (merge2,      out_lowres,  [('out', 'out_signal')]),
        (msk1,        out_lowres,  [('out_file', 'out_mask')]),
        (tpmmsk1,     out_lowres,  [('out1', 'out_tpms')]),
        (dist,        out_lowres,  [('outputnode.out_field', 'out_field'),
                                    ('outputnode.out_coeff', 'out_coeff'),
                                    ('outputnode.out_surfs', 'out_surfs')])
    ])
    return wf


def _bin_n_msk(in_files):
    import nibabel as nb
    import numpy as np
    import os.path as op

    in_files = np.atleast_1d(in_files).tolist()
    nii = [nb.load(f) for f in in_files]
    data = np.array([im.get_data() for im in nii])
    data = np.sum(data, axis=0)
    msk = np.zeros_like(data, dtype=np.uint8)
    msk[data > 0.0] = 1
    hdr = nii[0].get_header().copy()
    hdr.set_data_dtype(np.uint8)
    hdr.set_data_shape(msk.shape)
    out_file = op.abspath('binarized.nii.gz')
    nb.Nifti1Image(msk, nii[0].get_affine(), hdr).to_filename(out_file)
    return out_file


def _half_voxsize(in_file):
    import nibabel as nb
    import numpy as np
    return tuple(np.array(nb.load(in_file).get_header().get_zooms()[:3]) * 2.0)
