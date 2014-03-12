#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: Oscar Esteban - code@oscaresteban.es
# @Date:   2014-03-12 13:20:04
# @Last Modified by:   Oscar Esteban
# @Last Modified time: 2014-03-12 16:28:39

import os
import os.path as op
from glob import glob
import warnings

import numpy as np
import nibabel as nib

from nipype.interfaces.base import (traits, TraitedSpec, CommandLine,
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
                               desc='size of control points grid', argstr="-g %s" )
    out_prefix = traits.Str( "def", desc='output files prefix', argstr="-o %s", usedefault=True,
                             mandatory=True )

class RandomBSplineDeformationOutputSpec( TraitedSpec ):
    out_file = InputMultiPath(File(exists=True))
    out_coeff = InputMultiPath(File(exists=True))
    out_field = InputMultiPath(File(exists=True))

class RandomBSplineDeformation( CommandLine ):
    """ Use ACWEReg bspline random deformation tool to generate
    a deformation field and apply the deformation to the target
    images.

    Example
    -------
    >>> bspline = RandomBSplineDeformation()
    >>> bspline.inputs.in_file = 'moving.nii'
    >>> bspline.inputs.grid_size = [ 5, 8, 6 ]
    >>> bspline.inputs.out_prefix = 'myprefix'
    >>> bspline.cmdline
    'bspline_field -g 5 8 6 -I moving.nii -o myprefix'
    """

    input_spec = RandomBSplineDeformationInputSpec
    output_spec = RandomBSplineDeformationOutputSpec
    _cmd = '/home/oesteban/workspace/ACWE-Reg/Debug/Applications/bspline_field'

    def _format_arg( self, name, spec, value ):
        if name == "grid_size":
            if isinstance( value, list ):
                if not len(value)==3 and not len(value)==1:
                    raise RuntimeError("length of grid-size should be one value or three")

        return super(RandomBSplineDeformation, self)._format_arg( name, spec, value )

    def _list_outputs( self ):
        out_prefix = self.inputs.out_prefix
        outputs = self.output_spec().get()

        outputs['out_file'] = [ op.abspath( '%s_resampled_%d.nii.gz' % ( out_prefix, i ))  for i in range(len(self.inputs.in_file)) ]
        outputs['out_coeff'] = [ op.abspath( '%s_coeffs_%d.nii.gz' % ( out_prefix, i ))  for i in range(3) ]
        outputs['out_field'] = [ op.abspath( '%s_field_cmp%d.nii.gz' % ( out_prefix, i ))  for i in range(3) ]

        return outputs

