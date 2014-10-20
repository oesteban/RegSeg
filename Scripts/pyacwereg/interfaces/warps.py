# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: Oscar Esteban - code@oscaresteban.es
# @Date:   2014-03-12 13:20:04
# @Last Modified by:   oesteban
# @Last Modified time: 2014-10-15 00:15:22

import os
import os.path as op
import numpy as np
import nibabel as nib

from nipype.interfaces.base import (traits, TraitedSpec, CommandLine,
                                    BaseInterface, CommandLineInputSpec,
                                    BaseInterfaceInputSpec, InputMultiPath,
                                    OutputMultiPath, File, isdefined,
                                    Undefined)
from nipype.interfaces.ants.base import (ANTSCommand, ANTSCommandInputSpec)

from nipype import logging
logger = logging.getLogger('interface')


class RandomBSplineDeformationInputSpec(ANTSCommandInputSpec):
    in_file = InputMultiPath(File(exists=True), mandatory=True, argstr="-I %s",
                             desc='image(s) to be deformed')
    in_coeff = InputMultiPath(File(exists=True), mandatory=False,
                              argstr="-C %s",
                              desc=('coefficient images (to re-use a '
                                    'deformation field)'))
    in_surfs = InputMultiPath(File(exists=True), mandatory=False,
                              argstr="-S %s",
                              desc='apply deformation field to surfaces')
    in_mask = File(exists=True, argstr='-M %s', desc='set a mask')
    apply_mask = traits.Bool(False, argstr='--mask-inputs',
                             desc=('True if inputs should be masked out using '
                                   'deformed mask after transformation'))
    grid_size_item_trait = traits.Int(10, usedefault=True)
    grid_size = traits.Either(grid_size_item_trait,
                              traits.List(grid_size_item_trait), default=10,
                              argstr="-g %s", xor=['in_coeff'],
                              usedefault=True,
                              desc='size of control points grid')
    out_prefix = traits.Str("def", argstr="-o %s", usedefault=True,
                            mandatory=True, desc='output files prefix')


class RandomBSplineDeformationOutputSpec(TraitedSpec):
    out_file = OutputMultiPath(File(exists=True, desc='warped input files'))
    out_coeff = OutputMultiPath(File(exists=True, desc='output coefficients'))
    out_field = File(exists=True, desc='output warping field')
    out_field_base = OutputMultiPath(
        File(exists=True, desc='output warping bspline field'))
    out_surfs = OutputMultiPath(File(desc='output warped surfaces'))
    out_mask = File(exists=True, desc='warped input mask')


class RandomBSplineDeformation(ANTSCommand):

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
    _cmd = 'bspline_field'

    def _format_arg(self, name, spec, value):
        if name == "grid_size":
            if isinstance(value, list):
                if not len(value) == 3 and not len(value) == 1:
                    raise RuntimeError(
                        "length of grid-size should be one value or three")

        return super(RandomBSplineDeformation,
                     self)._format_arg(name, spec, value)

    def _list_outputs(self):
        out_prefix = self.inputs.out_prefix
        outputs = self.output_spec().get()

        outputs['out_file'] = [op.abspath('%s_warped_%d.nii.gz' % (
            out_prefix, i)) for i in range(len(self.inputs.in_file))]
        outputs['out_field'] = op.abspath('%s_dispfield.nii.gz' % out_prefix)
        outputs['out_coeff'] = [
            op.abspath('%s_coeffs_%d.nii.gz' % (out_prefix, i))
            for i in range(3)]
        outputs['out_field_base'] = [
            op.abspath('%s_field_cmp%d.nii.gz' % (out_prefix, i))
            for i in range(3)]

        if isdefined(self.inputs.in_surfs):
            outputs['out_surfs'] = [op.abspath(
                '%s_warped_%d.vtk' % (out_prefix, i))
                for i in range(len(self.inputs.in_surfs))]

        if isdefined(self.inputs.in_mask):
            outputs['out_mask'] = op.abspath(
                '%s_mask_warped.nii.gz' % out_prefix)

        return outputs


class FieldBasedWarpInputSpec(ANTSCommandInputSpec):
    in_file = InputMultiPath(File(exists=True), mandatory=True, argstr="-I %s",
                             desc='image(s) to be deformed')
    in_field = File(exists=True, mandatory=False, argstr="-F %s",
                    desc='forward field', xor='in_inv_field')
    in_inv_field = File(exists=True, mandatory=False, argstr="-R %s",
                        desc='backward field', xor='in_field')
    in_mask = File(exists=True, argstr='-M %s', desc='set a mask')
    in_surf = InputMultiPath(File(exists=True), argstr="-S %s",
                             desc='surface(s) to be deformed')
    out_prefix = traits.Str("fbased", argstr="-o %s", usedefault=True,
                            mandatory=True, desc='output files prefix')
    grid_size_item_trait = traits.Int(10, usedefault=True)
    grid_size = traits.Either(grid_size_item_trait,
                              traits.List(grid_size_item_trait), default=10,
                              argstr="-g %s", xor=['in_coeff'],
                              usedefault=True,
                              desc='size of control points grid')


class FieldBasedWarpOutputSpec(TraitedSpec):
    out_file = OutputMultiPath(File(exists=True, desc='warped input files'))
    out_surf = OutputMultiPath(File(exists=True, desc='warped input surfaces'))
    out_mask = File(exists=True, desc='warped input mask')


class FieldBasedWarp(ANTSCommand):

    """ Use ACWEReg bspline random deformation tool to generate
    a deformation field and apply the deformation to the target
    images.

    Example
    -------
    >>> warp = FieldBasedWarp()
    >>> warp.inputs.in_file = 'moving.nii'
    >>> warp.inputs.in_field = 'field.nii'
    >>> warp.inputs.in_surf = [ 'lh.white.vtk', 'rh.white.vtk' ]
    >>> warp.inputs.out_prefix = 'myprefix'
    >>> warp.cmdline
    'warp_image -I moving.nii -F field.nii -S lh.white.vtk rh.white.vtk \
-o myprefix'
    """

    input_spec = FieldBasedWarpInputSpec
    output_spec = FieldBasedWarpOutputSpec
    _cmd = 'warp_image'

    def _list_outputs(self):
        out_prefix = self.inputs.out_prefix
        outputs = self.output_spec().get()

        outputs['out_file'] = [op.abspath('%s_warped_%d.nii.gz' % (
            out_prefix, i)) for i in range(len(self.inputs.in_file))]

        if isdefined(self.inputs.in_mask):
            outputs['out_mask'] = op.abspath(
                '%s_mask_warped.nii.gz' % out_prefix)

        if isdefined(self.inputs.in_surf):
            outputs['out_surf'] = [op.abspath(
                '%s_warped_%d.vtk' % (out_prefix, i))
                for i in range(len(self.inputs.in_surf))]
        return outputs


class InverseFieldInputSpec(BaseInterfaceInputSpec):
    in_field = File(exists=True, mandatory=True,
                    desc='Displacements field to be inverted')
    out_field = File(desc='Output file name')


class InverseFieldOutputSpec(TraitedSpec):
    out_field = File(exists=True,
                     desc='Inverted displacement field')


class InverseField(BaseInterface):

    """ Takes a displacements field as input and computes the
    inverse of each vector.
    """
    input_spec = InverseFieldInputSpec
    output_spec = InverseFieldOutputSpec
    _out_file = ''

    def _run_interface(self, runtime):
        import numpy as np
        import nibabel as nb
        import os.path as op

        im = nb.load(self.inputs.in_field)
        data = -1.0 * im.get_data()
        nii = nb.Nifti1Image(data, im.get_affine(), im.get_header())
        if not isdefined(self.inputs.out_field):
            fname, ext = op.splitext(op.basename(self.inputs.in_field))
            if ext == '.gz':
                fname, ext2 = op.splitext(fname)
                ext = ext2 + ext
            out_file = op.abspath(fname + '_inv' + ext)
        else:
            out_file = self.inputs.out_field

        self._out_file = out_file

        nb.save(nii, out_file)

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_field'] = self._out_file
        return outputs
