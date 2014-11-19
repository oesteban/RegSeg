#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2014-11-19 09:46:07
# @Last Modified by:   oesteban
# @Last Modified time: 2014-11-19 10:14:49
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

        axis = 'axial'
        if self.inputs.axis == 0:
            axis = 'sagittal'
        elif self.inputs.axis == 1:
            axis = 'coronal'
        out_path = os.getcwd()
        out_pattern = op.join(out_path, '%s*.png' % axis)
        outputs['out_files'] = sorted(glob(out_pattern))
        return outputs
