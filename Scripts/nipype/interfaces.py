#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: Oscar Esteban - code@oscaresteban.es
# @Date:   2014-03-12 13:20:04
# @Last Modified by:   Oscar Esteban
# @Last Modified time: 2014-03-12 13:45:37

import os
from glob import glob
import warnings

import numpy as np
import nibabel as nib

from nipype.interfaces.base import (traits, TraitedSpec, , CommandLine,
                                    CommandLineInputSpec, InputMultiPath, File,
                                    isdefined, Undefined )

from nipype.utils.filemanip import load_json, save_json, split_filename, fname_presuffix

warn = warnings.warn
warnings.filterwarnings('always', category=UserWarning)

class RandomBSplineDeformationInputSpec( CommandLineInputSpec ):
    in_file = InputMultiPath(File(exists=True), mandatory=True, argstr="-I %s",
              desc='image(s) to be deformed')
    in_coeff = InputMultiPath(File(exists=True), mandatory=False, argstr="-C %s",
               desc='coefficient images (to re-use a deformation field)')
    in_surfs = InputMultiPath(File(exists=True), mandatory=False, argstr="-S %s",
               desc='apply deformation field to surfaces')
    grid_size_item_trait = traits.Int(10, usedefault=True )
    grid_size = traits.Either( grid_size_item_trait, traits.List(grid_size_item_trait),
                               xor=['in_coeff'], default=10, usedefault=True,
                               desc='size of control points grid' )
    out_prefix =

class RandomBSplineDeformationOutputSpec( TratedSpec ):
    out_file =
    out_coeff =
    out_field =


class RandomBSplineDeformation( CommandLine ):
    """ Use ACWEReg bspline random deformation tool to generate
    a deformation field and apply the deformation to the target
    images.
    """
    input_spec = RandomBSplineDeformationInputSpec
    output_spec = RandomBSplineDeformationOutputSpec
    _cmd = '/home/oesteban/workspace/ACWE-Reg/Debug/Applications/bspline_field'

    def _list_outputs( self ):
        out_prefix = self.inputs.out_prefix

        outputs = self.output_spec().get()
        return outputs

