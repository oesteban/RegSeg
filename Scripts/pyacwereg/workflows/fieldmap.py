#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2015-01-15 15:00:48
# @Last Modified by:   oesteban
# @Last Modified time: 2015-01-15 15:40:22

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
    polyfit = pe.Node(PolyFit(degree=4), name='FitPolyOrder4')
    scale = pe.Node(fsl.maths.BinaryMaths(operation='mul'), name='PhaseScaled')

    # Setup ANTS and registration
    def _aslist(tname):
        import numpy as np
        return np.atleast_1d(tname).tolist()

    fmm2t1w = pe.Node(ants.Registration(output_warped_image=True),
                      name="FMm_to_T1w")
    fmm2t1w.inputs.transforms = ['Rigid'] * 2
    fmm2t1w.inputs.transform_parameters = [(1.0,)] * 2
    fmm2t1w.inputs.number_of_iterations = [[50], [20]]
    fmm2t1w.inputs.dimension = 3
    fmm2t1w.inputs.metric = ['Mattes', 'Mattes']
    fmm2t1w.inputs.metric_weight = [1.0] * 2
    fmm2t1w.inputs.radius_or_number_of_bins = [64, 64]
    fmm2t1w.inputs.sampling_strategy = ['Regular', 'Random']
    fmm2t1w.inputs.sampling_percentage = [None, 0.2]
    fmm2t1w.inputs.convergence_threshold = [1.e-5, 1.e-8]
    fmm2t1w.inputs.convergence_window_size = [20, 10]
    fmm2t1w.inputs.smoothing_sigmas = [[6.0], [2.0]]
    fmm2t1w.inputs.sigma_units = ['vox'] * 2
    fmm2t1w.inputs.shrink_factors = [[6], [1]]  # ,[1] ]
    fmm2t1w.inputs.use_estimate_learning_rate_once = [True] * 2
    fmm2t1w.inputs.use_histogram_matching = [True] * 2
    fmm2t1w.inputs.initial_moving_transform_com = 0
    fmm2t1w.inputs.collapse_output_transforms = True
    fmm2t1w.inputs.winsorize_upper_quantile = 0.995

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
        (bet,                dilate, [('mask_file', 'in_file')]),
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
        (inputnode,           scale, [('factor', 'operand_value')]),
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
