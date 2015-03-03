#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2015-03-03 15:25:44
# @Last Modified by:   oesteban
# @Last Modified time: 2015-03-03 15:29:39


def mrtrix_dti(name='MRTrix_DTI'):
    """
    A workflow for DTI reconstruction using the tensor fitting included with
    MRTrix.
    :inputs:
        * in_dwi: the input dMRI volume to be reconstructed
        * in_bvec: b-vectors file in FSL format
        * in_bval: b-values file in FSL format
        * in_mask: input whole-brain mask (dwi space)
    """
    from nipype.pipeline import engine as pe
    from nipype.interfaces import utility as niu
    from nipype.interfaces import mrtrix as mrt
    from nipype.interfaces import freesurfer as fs

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_bvec', 'in_bval', 'in_dwi', 'in_mask']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['fa', 'md']), name='outputnode')

    fsl2mrtrix = pe.Node(mrt.FSL2MRTrix(), name='fsl2mrtrix')
    dwi2tsr = pe.Node(mrt.DWI2Tensor(), name='DWI2Tensor')
    tsr2fa = pe.Node(mrt.Tensor2FractionalAnisotropy(), name='ComputeFA')
    fa2nii = pe.Node(mrt.MRConvert(extension='nii'), 'FA2Nifti')
    msk_fa = pe.Node(fs.ApplyMask(), name='MaskFA')
    tsr2adc = pe.Node(mrt.Tensor2ApparentDiffusion(), name='ComputeADC')
    adc2nii = pe.Node(mrt.MRConvert(extension='nii'), 'ADC2Nifti')
    msk_adc = pe.Node(fs.ApplyMask(), name='MaskADC')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,    fsl2mrtrix, [('in_bvec', 'bvec_file'),
                                    ('in_bval', 'bval_file')]),
        (inputnode,       dwi2tsr, [('in_dwi', 'in_file')]),
        (inputnode,        msk_fa, [('in_mask', 'mask_file')]),
        (inputnode,       msk_adc, [('in_mask', 'mask_file')]),
        (fsl2mrtrix,      dwi2tsr, [('encoding_file', 'encoding_file')]),
        (dwi2tsr,          tsr2fa, [('tensor', 'in_file')]),
        (tsr2fa,           fa2nii, [('FA', 'in_file')]),
        (fa2nii,           msk_fa, [('converted', 'in_file')]),
        (dwi2tsr,         tsr2adc, [('tensor', 'in_file')]),
        (tsr2adc,         adc2nii, [('ADC', 'in_file')]),
        (adc2nii,         msk_adc, [('converted', 'in_file')]),
        (msk_fa,       outputnode, [('out_file', 'fa')]),
        (msk_adc,      outputnode, [('out_file', 'md')])
    ])

    return wf
