#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2015-01-15 15:00:48
# @Last Modified by:   oesteban
# @Last Modified time: 2015-02-03 16:43:19

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces import freesurfer as fs
from nipype.interfaces import fsl
from nipype.interfaces.io import JSONFileGrabber, JSONFileSink
from nipype.interfaces import ants
from nipype.algorithms.misc import AddNoise, PolyFit

from pysdcev.interfaces.misc import SigmoidFilter
import pysdcev.utils as pu


def bmap_registration(name="Bmap_Registration"):
    """
    A workflow to register a source B0 map to the T1w image of a real subject.
    """
    workflow = pe.Workflow(name=name)

    # Setup i/o
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['mag', 'pha', 't1w_brain', 'dwi_mask', 'poly_mask',
                'factor']), name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['magnitude', 'wrapped', 'unwrapped', 'mag_brain']),
        name='outputnode')

    # Setup initial nodes
    fslroi = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='GetFirst')

    mag2RAS = pe.Node(fs.MRIConvert(out_type="niigz", out_orientation="RAS"),
                      name='MagToRAS')
    pha2RAS = pe.Node(fs.MRIConvert(out_type="niigz", out_orientation="RAS"),
                      name='PhaToRAS')

    n4 = pe.Node(ants.N4BiasFieldCorrection(dimension=3), name='Bias')
    bet = pe.Node(fsl.BET(frac=0.4, mask=True), name='BrainExtraction')

    enh_mag = pe.Node(SigmoidFilter(upper_perc=98.8, lower_perc=10.0),
                      name='enh_mag')
    enh_t1w = pe.Node(SigmoidFilter(upper_perc=78.0, lower_perc=15.0),
                      name='enh_t1w')

    pha2rads = pe.Node(niu.Function(input_names=['in_file'],
                                    output_names=['out_file'],
                                    function=to_rads), name='Phase2rads')
    smooth = pe.Node(niu.Function(
        function=denoise, input_names=['in_file', 'in_mask'],
        output_names=['out_file']), name='PhaseDenoise')

    prelude = pe.Node(fsl.PRELUDE(process3d=True), name='PhaseUnwrap')

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

    polyfit = pe.Node(PolyFit(degree=2, padding=5), name='FitPoly')
    scale = pe.Node(niu.Function(
        function=scale_range, output_names=['out_file'],
        input_names=['in_file', 'value', 'in_mask']), name='ScaleBmap')
    scale.inputs.value = 1.8
    # scale = pe.Node(niu.Function(
    #     function=scale_like, output_names=['out_file'],
    #     input_names=['in_file', 'reference', 'in_mask']), name='ScaleBmap')

    median = pe.Node(niu.Function(
        input_names=['in_file'], output_names=['out_file'],
        function=pu.median_filter), name='PhaseMedian')
    demean = pe.Node(niu.Function(
        input_names=['in_file', 'in_mask'], output_names=['out_file'],
        function=pu.demean), name='PhaseDemean')
    addnoise = pe.Node(AddNoise(snr=30), name='PhaseAddNoise')
    wrap_pha = pe.Node(niu.Function(
        input_names=['in_file'], output_names=['out_file'],
        function=rads_ph_wrap), name='PhaseWrap')

    mwrapped = pe.Node(niu.Merge(2), name='MergeWrapped')
    munwrapped = pe.Node(niu.Merge(2), name='MergeUnwrapped')

    workflow.connect([
        (inputnode,        binarize, [('t1w_brain', 'in_file')]),
        (inputnode,          fslroi, [('mag', 'in_file')]),
        (inputnode,         pha2RAS, [('pha', 'in_file')]),
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

        # Transforms
        (inputnode,       warpPhase, [('t1w_brain', 'reference_image')]),
        (pha2RAS,          pha2rads, [('out_file', 'in_file')]),
        (bet,                smooth, [('mask_file', 'in_mask')]),
        (pha2rads,           smooth, [('out_file', 'in_file')]),
        (smooth,            prelude, [('out_file', 'phase_file')]),
        (n4,                prelude, [('output_image', 'magnitude_file')]),
        (prelude,         warpPhase, [
         ('unwrapped_phase_file', 'input_image')]),
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
        (regrid_pha,        polyfit, [('out_file', 'in_file')]),
        (inputnode,         polyfit, [('poly_mask', 'in_mask')]),
        (polyfit,             scale, [('out_file', 'in_file')]),
        # (regrid_pha,          scale, [('out_file', 'reference')]),
        (inputnode,           scale, [('poly_mask', 'in_mask')]),
        (scale,              median, [('out_file', 'in_file')]),
        (median,             demean, [('out_file', 'in_file')]),
        (inputnode,          demean, [('dwi_mask', 'in_mask')]),
        (demean,           addnoise, [('out_file', 'in_file')]),
        (inputnode,        addnoise, [('dwi_mask', 'in_mask')]),
        (addnoise,         wrap_pha, [('out_file', 'in_file')]),
        (regrid_bmg,     munwrapped, [('out_file', 'in1')]),
        (demean,         munwrapped, [('out_file', 'in2')]),
        (regrid_bmg,       mwrapped, [('out_file', 'in1')]),
        (wrap_pha,         mwrapped, [('out_file', 'in2')]),
        (regrid_mag,     outputnode, [('out_file', 'magnitude')]),
        (regrid_bmg,     outputnode, [('out_file', 'mag_brain')]),
        (munwrapped,     outputnode, [('out', 'unwrapped')]),
        (mwrapped,       outputnode, [('out', 'wrapped')]),
    ])
    return workflow


def denoise(in_file, in_mask, out_file=None):
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


def scale_range(in_file, value=2.0, in_mask=None, out_file=None):
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
    idata -= (factor * ip0) + rp0

    nb.Nifti1Image(
        idata, im.get_affine(), im.get_header()).to_filename(out_file)

    return out_file
