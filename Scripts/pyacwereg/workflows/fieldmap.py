#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2015-01-15 15:00:48
# @Last Modified by:   oesteban
# @Last Modified time: 2015-03-11 10:24:25

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces import freesurfer as fs
from nipype.interfaces import fsl
from nipype.interfaces.io import JSONFileGrabber, JSONFileSink
from nipype.interfaces import ants
from nipype.algorithms.misc import AddNoise
from nipype.workflows.dmri.fsl import utils as nwu

from pyacwereg.interfaces.dmri import PhaseUnwrap
from pyacwereg.interfaces.utility import SigmoidFilter


def bmap_registration(name="Bmap_Registration"):
    """
    A workflow to register a source B0 map to the T1w image of a real subject.
    """
    workflow = pe.Workflow(name=name)

    # Setup i/o
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['mag', 'pha', 't1w_brain', 'dwi_mask']),
        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['magnitude', 'wrapped', 'unwrapped', 'mag_brain']),
        name='outputnode')

    # Setup initial nodes
    fslroi = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='GetFirst')

    mag2RAS = pe.Node(fs.MRIConvert(out_type="niigz", out_orientation="RAS"),
                      name='MagToRAS')
    unwrap = pe.Node(PhaseUnwrap(), name='PhaseUnwrap')
    pha2RAS = pe.Node(fs.MRIConvert(out_type="niigz", out_orientation="RAS"),
                      name='PhaToRAS')

    n4 = pe.Node(ants.N4BiasFieldCorrection(dimension=3), name='Bias')
    bet = pe.Node(fsl.BET(frac=0.4, mask=True), name='BrainExtraction')

    enh_mag = pe.Node(SigmoidFilter(upper_perc=98.8, lower_perc=10.0),
                      name='enh_mag')
    enh_t1w = pe.Node(SigmoidFilter(upper_perc=78.0, lower_perc=15.0),
                      name='enh_t1w')

    # Setup ANTS and registration
    def _aslist(tname):
        import numpy as np
        return np.atleast_1d(tname).tolist()

    fmm2t1w = pe.Node(ants.Registration(output_warped_image=True),
                      name="FMm_to_T1w")
    fmm2t1w.inputs.transforms = ['Rigid'] * 2
    fmm2t1w.inputs.transform_parameters = [(1.0,)] * 2
    fmm2t1w.inputs.number_of_iterations = [[250], [100]]
    fmm2t1w.inputs.dimension = 3
    fmm2t1w.inputs.metric = ['Mattes', 'Mattes']
    fmm2t1w.inputs.metric_weight = [1.0] * 2
    fmm2t1w.inputs.radius_or_number_of_bins = [64, 64]
    fmm2t1w.inputs.sampling_strategy = ['Regular', 'Random']
    fmm2t1w.inputs.sampling_percentage = [None, 0.1]
    fmm2t1w.inputs.convergence_threshold = [1.e-5, 1.e-7]
    fmm2t1w.inputs.convergence_window_size = [10, 5]
    fmm2t1w.inputs.smoothing_sigmas = [[6.0], [2.0]]
    fmm2t1w.inputs.sigma_units = ['vox'] * 2
    fmm2t1w.inputs.shrink_factors = [[6], [1]]  # ,[1] ]
    fmm2t1w.inputs.use_estimate_learning_rate_once = [True] * 2
    fmm2t1w.inputs.use_histogram_matching = [True] * 2
    fmm2t1w.inputs.initial_moving_transform_com = 0
    fmm2t1w.inputs.collapse_output_transforms = True

    binarize = pe.Node(fs.Binarize(min=0.1), name='BinT1')

    warpPhase = pe.Node(ants.ApplyTransforms(
        dimension=3, interpolation='BSpline'), name='WarpPhase')
    warpMag = pe.Node(ants.ApplyTransforms(
        dimension=3, interpolation='BSpline'), name='WarpMag')

    # Final regrids and phase re-wrapping
    regrid_mag = pe.Node(fs.MRIConvert(
        resample_type='cubic', out_datatype='float'), name='Regrid_mag')
    regrid_bmg = pe.Node(fs.MRIConvert(
        resample_type='cubic', out_datatype='float'), name='Regrid_mag_brain')
    regrid_pha = pe.Node(fs.MRIConvert(
        resample_type='cubic', out_datatype='float'), name='Regrid_pha')

    # denoise = pe.Node(niu.Function(
    #     input_names=['in_file', 'in_mask'], output_names=['out_file'],
    #     function=filter_fmap), name='SmoothBmap')

    addnoise = pe.Node(AddNoise(snr=30), name='PhaseAddNoise')
    wrap_pha = pe.Node(niu.Function(
        input_names=['in_file'], output_names=['out_file'],
        function=rads_ph_wrap), name='PhaseWrap')

    mwrapped = pe.Node(niu.Merge(2), name='MergeWrapped')
    munwrapped = pe.Node(niu.Merge(2), name='MergeUnwrapped')

    workflow.connect([
        (inputnode,        binarize, [('t1w_brain', 'in_file')]),
        (inputnode,          fslroi, [('mag', 'in_file')]),
        # Connect first nodes
        (fslroi,            mag2RAS, [('roi_file', 'in_file')]),
        (mag2RAS,                n4, [('out_file', 'input_image')]),
        (n4,                    bet, [('output_image', 'in_file')]),
        (inputnode,         enh_t1w, [('t1w_brain', 'in_file')]),
        (binarize,          enh_t1w, [('binary_file', 'in_mask')]),
        (n4,                enh_mag, [('output_image', 'in_file')]),
        (bet,               enh_mag, [('mask_file', 'in_mask')]),
        # ANTs
        (enh_t1w,           fmm2t1w, [('out_file', 'fixed_image')]),
        (enh_mag,           fmm2t1w, [('out_file', 'moving_image')]),
        (binarize,          fmm2t1w, [('binary_file', 'fixed_image_mask')]),
        (bet,               fmm2t1w, [('mask_file', 'moving_image_mask')]),

        # Unwrap
        (inputnode,          unwrap, [('pha', 'in_file')]),
        (unwrap,            pha2RAS, [('out_file', 'in_file')]),

        # Transforms
        (inputnode,       warpPhase, [('t1w_brain', 'reference_image')]),
        (pha2RAS,         warpPhase, [('out_file', 'input_image')]),
        (fmm2t1w,    warpPhase, [
            ('forward_transforms', 'transforms'),
            ('forward_invert_flags', 'invert_transform_flags')]),
        (inputnode,         warpMag, [('t1w_brain', 'reference_image')]),
        (n4,                warpMag, [('output_image', 'input_image')]),
        (fmm2t1w,           warpMag, [
            ('forward_transforms', 'transforms'),
            ('forward_invert_flags', 'invert_transform_flags')]),
        (warpMag,        regrid_mag, [('output_image', 'in_file')]),
        (inputnode,      regrid_mag, [('dwi_mask', 'reslice_like')]),
        (fmm2t1w,        regrid_bmg, [('warped_image', 'in_file')]),
        (inputnode,      regrid_bmg, [('dwi_mask', 'reslice_like')]),
        (warpPhase,      regrid_pha, [('output_image', 'in_file')]),
        (inputnode,      regrid_pha, [('dwi_mask', 'reslice_like')]),
        (regrid_pha,        addnoise, [('out_file', 'in_file')]),
        # (regrid_pha,        denoise, [('out_file', 'in_file')]),
        # (inputnode,         denoise, [('dwi_mask', 'in_mask')]),
        # (denoise,          addnoise, [('out_file', 'in_file')]),
        (inputnode,        addnoise, [('dwi_mask', 'in_mask')]),
        (addnoise,         wrap_pha, [('out_file', 'in_file')]),
        (regrid_bmg,     munwrapped, [('out_file', 'in1')]),
        (regrid_pha,     munwrapped, [('out_file', 'in2')]),
        # (denoise,        munwrapped, [('out_file', 'in2')]),
        (regrid_bmg,       mwrapped, [('out_file', 'in1')]),
        (wrap_pha,         mwrapped, [('out_file', 'in2')]),
        (regrid_mag,     outputnode, [('out_file', 'magnitude')]),
        (regrid_bmg,     outputnode, [('out_file', 'mag_brain')]),
        (munwrapped,     outputnode, [('out', 'unwrapped')]),
        (mwrapped,       outputnode, [('out', 'wrapped')]),
    ])
    return workflow


def vsm_fmb(name='VSM_FMB', phase_unwrap=True, fmb_params={}):
    """
    Workflow that uses `FUGUE <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE>`_
    to compute the voxel-shift map (VSM) from the corresponding B-fieldmap in
    EPI data.
    """
    wf = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_bmap', 'in_mask', 'echospacing', 'delta_te', 'acc_factor',
                'enc_dir']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['vsm', 'dfm', 'dfm_inv', 'jac', 'jac_inv']), name='outputnode')

    selmag = pe.Node(niu.Select(index=0), name='SelectMag')
    selpha = pe.Node(niu.Select(index=1), name='SelectPha')

    magmsk = pe.Node(fs.ApplyMask(), name='mask_mag')
    eff_es = pe.Node(niu.Function(
        input_names=['echospacing', 'acc_factor'], function=_eff_es_seg,
        output_names=['echospacing']), name='GetEffEcho')
    phcache = pe.Node(niu.IdentityInterface(fields=['pha_unwrapped']),
                      name='PhCache')

    rad2rsec = pe.Node(niu.Function(
        input_names=['in_file', 'delta_te'], output_names=['out_file'],
        function=nwu.rads2radsec), name='ToRadSec')
    pre_fugue = pe.Node(fsl.FUGUE(save_unmasked_fmap=True),
                        name='PreliminaryFugue')
    demean = pe.Node(niu.Function(
        function=nwu.demean_image, input_names=['in_file', 'in_mask'],
        output_names=['out_file']), name='DemeanFmap')
    cleanup = nwu.cleanup_edge_pipeline()
    addvol = pe.Node(niu.Function(
        function=nwu.add_empty_vol, input_names=['in_file'],
        output_names=['out_file']), name='AddEmptyVol')
    vsm = pe.Node(fsl.FUGUE(save_unmasked_shift=True, **fmb_params),
                  name="Compute_VSM")

    dfm = process_vsm()
    dfm.inputs.inputnode.scaling = 1.0

    wf.connect([
        (inputnode,   selmag,     [('in_bmap', 'inlist')]),
        (inputnode,   selpha,     [('in_bmap', 'inlist')]),
        (inputnode,   magmsk,     [('in_mask', 'mask_file')]),
        (inputnode,   eff_es,     [('echospacing', 'echospacing'),
                                   ('acc_factor', 'acc_factor')]),
        (selmag,      magmsk,     [('out', 'in_file')]),
        (phcache,     rad2rsec,   [('pha_unwrapped', 'in_file')]),
        (inputnode,   rad2rsec,   [('delta_te', 'delta_te')]),
        (rad2rsec,    pre_fugue,  [('out_file', 'fmap_in_file')]),
        (inputnode,   pre_fugue,  [('in_mask', 'mask_file')]),
        (pre_fugue,   demean,     [('fmap_out_file', 'in_file')]),
        (inputnode,   demean,     [('in_mask', 'in_mask')]),
        (demean,      cleanup,    [('out_file', 'inputnode.in_file')]),
        (inputnode,   cleanup,    [('in_mask', 'inputnode.in_mask')]),
        (cleanup,     addvol,     [('outputnode.out_file', 'in_file')]),
        (addvol,      vsm,        [('out_file', 'fmap_in_file')]),
        (eff_es,      vsm,        [('echospacing', 'dwell_time')]),
        (inputnode,   vsm,        [('in_mask', 'mask_file'),
                                   ('delta_te', 'asym_se_time')]),
        (vsm,         outputnode, [('shift_out_file', 'vsm')]),
        (vsm,         dfm,        [('shift_out_file', 'inputnode.vsm')]),
        (inputnode,   dfm,        [('in_mask', 'inputnode.reference'),
                                   ('enc_dir', 'inputnode.enc_dir')]),
        (dfm,         outputnode, [('outputnode.dfm', 'dfm'),
                                   ('outputnode.jacobian', 'jac'),
                                   ('outputnode.inv_dfm', 'dfm_inv'),
                                   ('outputnode.inv_jacobian', 'jac_inv')])
    ])

    if phase_unwrap:
        prelude = pe.Node(fsl.PRELUDE(process3d=True), name='PhaseUnwrap')
        wf.connect([
            (selmag,   prelude, [('out', 'magnitude_file')]),
            (selpha,   prelude, [('out', 'phase_file')]),
            (magmsk,   prelude, [('out_file', 'mask_file')]),
            (prelude,  phcache, [('unwrapped_phase_file', 'pha_unwrapped')])
        ])
    else:
        wf.connect(selpha, 'out', phcache, 'pha_unwrapped')
    return wf


def process_vsm(name='VSM2DFM'):
    """
    A workflow that computes the inverse and the determinant of jacobians
    map
    """
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['vsm', 'reference', 'scaling', 'enc_dir']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['dfm', 'inv_dfm', 'jacobian', 'inv_jacobian']),
        name='outputnode')

    vsm2dfm = nwu.vsm2warp()
    invwarp = pe.Node(fsl.InvWarp(relative=True), name='InvWarp')
    coeffs = pe.Node(fsl.WarpUtils(out_format='spline'),
                     name='CoeffComp')
    jacobian = pe.Node(fsl.WarpUtils(write_jacobian=True),
                       name='JacobianComp')
    invcoef = pe.Node(fsl.WarpUtils(out_format='spline'),
                      name='InvCoeffComp')
    invjac = pe.Node(fsl.WarpUtils(write_jacobian=True),
                     name='InvJacobianComp')
    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, vsm2dfm,    [('reference', 'inputnode.in_ref'),
                                 ('enc_dir', 'inputnode.enc_dir'),
                                 ('vsm', 'inputnode.in_vsm'),
                                 ('scaling', 'inputnode.scaling')]),
        (vsm2dfm,   invwarp,    [('outputnode.out_warp', 'warp')]),
        (inputnode, invwarp,    [('reference', 'reference')]),
        (vsm2dfm,   coeffs,     [('outputnode.out_warp', 'in_file')]),
        (inputnode, coeffs,     [('reference', 'reference')]),
        (inputnode, jacobian,   [('reference', 'reference')]),
        (coeffs,    jacobian,   [('out_file', 'in_file')]),
        (vsm2dfm,   outputnode, [('outputnode.out_warp', 'dfm')]),
        (invwarp,   outputnode, [('inverse_warp', 'inv_dfm')]),
        (jacobian,  outputnode, [('out_jacobian', 'jacobian')]),
        (invwarp,   invcoef,    [('inverse_warp', 'in_file')]),
        (inputnode, invcoef,    [('reference', 'reference')]),
        (inputnode, invjac,     [('reference', 'reference')]),
        (invcoef,   invjac,     [('out_file', 'in_file')]),
        (invjac,    outputnode, [('out_jacobian', 'inv_jacobian')])
    ])
    return wf


def _eff_es_seg(echospacing, acc_factor):
    return float(echospacing / (1.0 * acc_factor))


def median_f(in_file, in_mask, out_file=None):
    import nibabel as nb
    import os.path as op
    import numpy as np
    from scipy import ndimage as nd

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_smoothed.nii.gz' % fname)

    im = nb.load(in_file)
    noisy = im.get_data().astype(np.float32)
    msk = nb.load(in_mask).get_data()
    data = nd.median_filter(noisy, 5)
    data[msk > 0.] = noisy[msk > 0.]
    nb.Nifti1Image(data, im.get_affine(), im.get_header()).to_filename(
        out_file)
    return out_file


def to_rads(in_file, out_file=None):
    import nibabel as nb
    import os.path as op
    import numpy as np
    import math

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_rads.nii.gz' % fname)

    im = nb.load(in_file)
    hdr = im.get_header().copy()
    data = im.get_data().astype(np.float32)
    imin = data.min()
    imax = data.max()

    data = (2.0 * math.pi * (data - imin) / (imax - imin)) - math.pi
    hdr.set_data_dtype(np.float32)
    hdr['datatype'] = 16
    nb.Nifti1Image(data, im.get_affine(), hdr).to_filename(out_file)
    return out_file


def to_rad_sec(in_file, mask_file=None, demean=False, delta_te=2.46e-3,
               out_file=None):
    import nibabel as nb
    import os.path as op
    import numpy as np
    import math

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_radsec.nii.gz' % fname)

    im = nb.load(in_file)
    data = im.get_data().astype(np.float32) * (1.0 / delta_te)

    if demean:
        if mask_file is None:
            raise RuntimeError('A mask should be supplied to demean')

        msk = nb.load(mask_file).get_data()
        msk[msk <= 0] = 0
        msk[msk > 0] = 1

        data = np.ma.array(data, mask=1 - msk)
        mval = np.ma.median(data)
        data = data - mval

    nb.Nifti1Image(data, im.get_affine(),
                   im.get_header()).to_filename(out_file)
    return out_file


def rads_ph_wrap(in_file, out_file=None, orange=None):
    import nibabel as nb
    import os.path as op
    import numpy as np
    import math

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_wrapped.nii.gz' % fname)

    im = nb.load(in_file)
    imdata = im.get_data().astype(np.float32) + math.pi
    periods = 2.0 * math.pi * (imdata // (2.0 * math.pi))
    imdata = imdata - periods
    hdr = im.get_header().copy()

    if orange is not None:
        imdata = (orange * (imdata / (2.0 * math.pi)) -
                  (0.5 * orange - 1)).astype(np.int16)
        hdr.set_data_dtype(np.int16)
        hdr['datatype'] = 4
    else:
        imdata = imdata - math.pi

    nb.Nifti1Image(imdata, im.get_affine(),
                   hdr).to_filename(out_file)
    return out_file


def siemens_ph_wrap(in_file, out_file=None, irange=8192, orange=4096.00):
    import nibabel as nb
    import os.path as op
    import numpy as np
    import math

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_wrapped.nii.gz' % fname)

    im = nb.load(in_file)
    imdata = im.get_data().astype(np.float32)
    imdata = (imdata - np.median(imdata)).astype(np.int16)

    imax = int(irange * 0.5)
    imin = int(irange * -0.5) + 1

    imdata[imdata > imax] = imdata[imdata > imax] - irange
    imdata[imdata < imin] = imdata[imdata < imin] + irange

    imdata = imdata.astype(np.float32) / (1.0 * irange)
    imdata = imdata * orange

    hdr = im.get_header().copy()
    hdr.set_data_dtype(np.int16)
    hdr['datatype'] = 4

    nii = nb.Nifti1Image(imdata.astype(np.int16), im.get_affine(), hdr)
    nb.save(nii, out_file)
    return out_file


def _b0_field(b0_axis, b0_strength):
    import numpy as np
    b0_dir = np.array([0.0, 0.0, 1.0]) * float(b0_strength)
    if b0_axis != 'z':
        b0_dir = np.roll(b0_dir, 1)
        if b0_axis == 'y':
            b0_dir = np.roll(b0_dir, 1)

    return tuple(b0_dir.tolist())


def bmap2vsm(in_file, echospacing, acc_factor, enc_dir='y', in_mask=None,
             out_file=None):
    import numpy as np
    import nibabel as nb
    import os.path as op
    import math

    im = nb.load(in_file)

    if im.get_data().ndim == 4:
        im = nb.funcs.four_to_three(im)[0]

    sizes = im.get_shape()

    N = sizes[1]

    if enc_dir[0] == 'x':
        N = sizes[0]
    elif enc_dir[0] == 'z':
        N = sizes[2]

    direction = -1.0

    if len(enc_dir) == 2 and enc_dir[1] == '-':
        direction = 1.0

    data = im.get_data()

    if in_mask is not None:
        msk = nb.load(in_mask).get_data()
        mskdata = np.ma.array(data, mask=1 - msk)
        mval = np.ma.median(mskdata)
        data = data - mval

    eff = float(echospacing / (1.0 * acc_factor))
    data = direction * 42.57748e6 * eff * N * data

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_vsm.nii.gz' % fname)

    nb.Nifti1Image(data.astype(np.float32), im.get_affine(),
                   im.get_header()).to_filename(out_file)

    return out_file


def bmap2phasediff(in_file, delta_te, in_mask=None, out_file=None):
    import numpy as np
    import nibabel as nb
    import os.path as op
    import math

    im = nb.load(in_file)
    data = im.get_data()

    if in_mask is not None:
        msk = nb.load(in_mask).get_data()
        mskdata = np.ma.array(data, mask=1 - msk)
        mval = np.ma.median(mskdata)
        data = data - mval

    data = (2.0 * math.pi) * data * 42.57748e6 * float(delta_te)

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_phasemap.nii.gz' % fname)

    nb.Nifti1Image(data.astype(np.float32), im.get_affine(),
                   im.get_header()).to_filename(out_file)

    return out_file


def phasediff2siemens(in_file, out_file=None):
    import numpy as np
    import nibabel as nb
    import math
    import os.path as op

    im = nb.load(in_file)

    data = (2048 / (2.0 * math.pi) * im.get_data()).astype(np.int16)
    data[data > 2048] = data[data > 2048] - 4096
    data[data < -2047] = data[data < -2047] + 4096

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_siemens.nii.gz' % fname)

    hdr = im.get_header().copy()
    hdr.set_data_dtype(np.int16)
    im1 = nb.Nifti1Image(data.astype(np.int16), im.get_affine(), hdr)
    im2 = nb.Nifti1Image(np.zeros_like(data), im.get_affine(), hdr)
    nb.funcs.concat_images([im1, im2]).to_filename(out_file)

    return out_file


def scale_like(in_file, reference, in_mask=None, out_file=None):
    import numpy as np
    import nibabel as nb
    import os.path as op

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_slike.nii.gz' % fname)

    im = nb.load(in_file)
    idata = im.get_data()
    rdata = nb.load(reference).get_data()

    mask = np.ones_like(idata)
    if in_mask is not None:
        mask = nb.load(in_mask).get_data()
        mask[mask > 0.0] = 1
        mask[mask < 0.0] = 0

    rp0 = np.percentile(rdata[mask > 0.0], 0.5)
    rp1 = np.percentile(rdata[mask > 0.0], 99.5)
    ip0 = np.percentile(idata[mask > 0.0], 0.5)
    ip1 = np.percentile(idata[mask > 0.0], 99.5)

    factor = (rp1 - rp0) / (ip1 - ip0)
    idata *= factor
    idata -= (factor * ip0) + rp0

    nb.Nifti1Image(
        idata, im.get_affine(), im.get_header()).to_filename(out_file)

    return out_file


def scale_range(in_file, value=1.0, in_mask=None, out_file=None):
    import numpy as np
    import nibabel as nb
    import os.path as op

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_range.nii.gz' % fname)

    im = nb.load(in_file)
    idata = im.get_data()

    mask = np.ones_like(idata)
    if in_mask is not None:
        mask = nb.load(in_mask).get_data()
        mask[mask > 0.0] = 1
        mask[mask < 0.0] = 0

    rp0 = -1.0 * value
    rp1 = value
    ip0 = np.percentile(idata[mask > 0.0], 1.5)
    ip1 = np.percentile(idata[mask > 0.0], 98.5)

    factor = (rp1 - rp0) / (ip1 - ip0)
    idata *= factor
    idata -= np.median(idata[mask > 0.0])

    nb.Nifti1Image(
        idata, im.get_affine(), im.get_header()).to_filename(out_file)

    return out_file


def filter_fmap(in_file, in_mask=None, out_file=None):
    from pyacwereg.filters import wavelets_denoise, laplacian_filter
    import numpy as np
    import nibabel as nb
    import os.path as op
    from math import pi
    from scipy.ndimage import median_filter

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_filtered.nii.gz' % fname)

    filtered = wavelets_denoise(laplacian_filter(in_file, in_mask))

    im = nb.load(filtered)
    result = im.get_data()
    result = median_filter(result, 10)
    result -= np.median(result)
    result *= (pi / np.percentile(result, 99.97))
    nb.Nifti1Image(result, im.get_affine(),
                   im.get_header()).to_filename(out_file)

    return out_file
