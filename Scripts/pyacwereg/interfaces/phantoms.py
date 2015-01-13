#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2014-10-23 14:45:06
# @Last Modified by:   oesteban
# @Last Modified time: 2015-01-13 13:01:17
import os
import os.path as op
import nibabel as nb
import numpy as np

from nipype.interfaces.base import (BaseInterface, traits, TraitedSpec, File,
                                    InputMultiPath, OutputMultiPath,
                                    BaseInterfaceInputSpec, isdefined,
                                    DynamicTraitedSpec, Directory,
                                    CommandLine, CommandLineInputSpec)

from pyacwereg import misc as pm
from nipype import logging
iflogger = logging.getLogger('interface')


class PhantomInputSpec(BaseInterfaceInputSpec):
    shape = traits.Enum('gyrus', 'box', 'L', 'ball', usedefault=True,
                        desc='Phantom shape')
    cortex = traits.Bool(True, usedefault=True,
                         desc='Generate a crust mimicking cortical GM')
    out_file = File('phantom_model.nii.gz', usedefault=True,
                    desc='output file name')
    matrix_size_item_trait = traits.Int(101, usedefault=True)
    matrix_size = traits.Either(matrix_size_item_trait,
                                traits.List(matrix_size_item_trait),
                                default=101, usedefault=True,
                                desc='image matrix')


class PhantomOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='output file name')
    out_mask = File(exists=True, desc='output file name')


class Phantom(BaseInterface):

    """
    Returns a phantom model
    """
    input_spec = PhantomInputSpec
    output_spec = PhantomOutputSpec

    def _run_interface(self, runtime):
        size = np.atleast_1d(self.inputs.matrix_size).tolist()
        if len(size) == 1:
            size = [size[0]] * 3

        data = pm.genShape(self.inputs.shape, datashape=size,
                           cortex=self.inputs.cortex)
        nii = pm.genNiftiVol(data)
        nii.to_filename(op.abspath(self.inputs.out_file))

        mask = (np.ones_like(data[0]) - data[0]).astype(np.uint8)
        hdr = nii.get_header().copy()
        hdr.set_data_dtype(np.uint8)
        hdr.set_data_shape(mask.shape)
        hdr['data_type'] = 2
        nb.Nifti1Image(mask, nii.get_affine(), hdr).to_filename(
            self._get_outmsk())
        return runtime

    def _get_outmsk(self):
        out_mask, ext = op.splitext(op.abspath(self.inputs.out_file))
        if ext == '.gz':
            out_mask, ext2 = op.splitext(out_mask)
            ext = ext2 + ext
        out_mask = out_mask + '_mask' + ext
        return out_mask

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = op.abspath(self.inputs.out_file)
        outputs['out_mask'] = self._get_outmsk()
        return outputs


class SimulateSMRIInputSpec(CommandLineInputSpec):
    frac_csf = File(exists=True, argstr='--csf_vf %s',
                    desc='CSF volume fraction')
    frac_wm = File(exists=True, mandatory=True, argstr='--wm_vf %s',
                   desc='WM volume fraction')
    frac_gm = File(exists=True, mandatory=True, argstr='--gm_vf %s',
                   desc='GM volume fraction')
    out_dir = Directory(argstr='-o %s', desc='output prefix')
    snr = traits.Float(0.0, usedefault=True, argstr='--snr %f',
                       desc='SNR of output images')


class SimulateSMRIOutputSpec(TraitedSpec):
    out_t1w = File(exists=True, desc='output file name')
    out_t2w = File(exists=True, desc='output file name')


class SimulateSMRI(CommandLine):

    """
    Returns a phantom model
    """
    input_spec = SimulateSMRIInputSpec
    output_spec = SimulateSMRIOutputSpec
    cmd = 'phantomas_struct_fiberless'

    def _list_outputs(self):
        outputs = self._outputs().get()

        if isdefined(self.inputs.out_dir):
            out_dir = op.abspath(self.inputs.out_dir)
        else:
            out_dir = os.getcwd()

        outputs['out_t1w'] = op.join(out_dir, 't1_weighted.nii.gz')
        outputs['out_t2w'] = op.join(out_dir, 't2_weighted.nii.gz')
        return outputs


class DownsampleAveragingInputSpec(CommandLineInputSpec):
    in_file = File(exists=True, position=1, argstr='-i %s', mandatory=True,
                   desc='input binary file')
    out_file = File('downsampled.nii.gz', usedefault=True, argstr='-o %s')
    matrix_size_item_trait = traits.Int(51, usedefault=True)
    matrix_size = traits.Either(matrix_size_item_trait,
                                traits.Tuple(matrix_size_item_trait,
                                             matrix_size_item_trait,
                                             matrix_size_item_trait),
                                argstr='-g %s',
                                default=51, usedefault=True,
                                desc='image matrix')


class DownsampleAveragingOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='output file name')


class DownsampleAveraging(CommandLine):
    input_spec = DownsampleAveragingInputSpec
    output_spec = DownsampleAveragingOutputSpec
    cmd = 'regridavg'

    def _format_arg(self, name, spec, value):
        if name == 'matrix_size':
            size = np.atleast_1d(value).tolist()
            if len(size) == 1:
                size = [size[0]] * 3
            return spec.argstr % " ".join(['%d' % s for s in size])
        if name == 'out_file':
            return spec.argstr % op.abspath(value)
        return super(DownsampleAveraging, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = op.abspath(self.inputs.out_file)
        return outputs
