#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2014-11-19 09:46:07
# @Last Modified by:   oesteban
# @Last Modified time: 2015-01-13 12:35:51
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


class ExportSlicesInputSpec(CommandLineInputSpec):
    reference = File(exists=True, argstr='-i %s', mandatory=True,
                     desc=('reference image to show in background'))
    surfaces0 = InputMultiPath(
        File(exists=True), argstr='-S %s',
        desc=('vtk contours that will be overlaid on reference'))
    surfaces1 = InputMultiPath(
        File(exists=True), argstr='-R %s',
        desc=('vtk contours that will be overlaid on reference'))
    num_slices = traits.Int(14, argstr='-n %d', desc='total num. of slices')
    axis = traits.Enum(2, 0, 1, argstr='-a %d', usedefault=True,
                       desc='axis to cut through')
    all_axis = traits.Bool(argstr='-A', desc='slice through all axes')


class ExportSlicesOutputSpec(TraitedSpec):
    out_files = OutputMultiPath(File(exists=True), desc='output files')


class ExportSlices(CommandLine):

    """
    Export slices produces visualization through the selected axis
    of the input image and contours
    """
    input_spec = ExportSlicesInputSpec
    output_spec = ExportSlicesOutputSpec
    _cmd = 'slice_contours'

    def _list_outputs(self):
        from glob import glob
        outputs = self.output_spec().get()
        out_path = os.getcwd()
        out_pattern = op.join(out_path, '*.png')
        outputs['out_files'] = sorted(glob(out_pattern))
        return outputs


class Surf2VolInputSpec(CommandLineInputSpec):
    reference = File(exists=True, argstr='-R %s', mandatory=True,
                     desc=('reference image grid'))
    surfaces = InputMultiPath(
        File(exists=True), argstr='-S %s', mandatory=True,
        desc=('vtk contours that will be mapped to volume'))
    out_prefix = traits.Str('surf2vol', argstr='-o %s', usedefault=True,
                            desc='output files prefix')


class Surf2VolOutputSpec(TraitedSpec):
    out_tpm = OutputMultiPath(File(exists=True),
                              desc='output tissue probability maps')
    out_seg = File(exists=True, desc='output segmentation')


class Surf2Vol(CommandLine):

    """
    Converts surface contours defining regions in space to volumes
    """
    input_spec = Surf2VolInputSpec
    output_spec = Surf2VolOutputSpec
    _cmd = 'binarize_surfaces'

    def _list_outputs(self):
        outputs = self.output_spec().get()
        out_path = op.abspath(self.inputs.out_prefix)

        outputs['out_tpm'] = [
            out_path + ('_tpm_cmp%d.nii.gz' % i)
            for i in range(len(self.inputs.surfaces))]
        outputs['out_seg'] = out_path + '_seg.nii.gz'
        return outputs
