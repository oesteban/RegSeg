# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: Oscar Esteban - code@oscaresteban.es
# @Date:   2014-03-12 13:20:04
# @Last Modified by:   oesteban
# @Last Modified time: 2014-10-22 12:30:33

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


class ACWERegInputGroupSpec(ANTSCommandInputSpec):
    # Functional options
    float_trait = traits.Either(None, traits.Float(1.0))
    int_trait = traits.Either(None, traits.Int(0))
    bool_trait = traits.Either(None, traits.Bool(False))

    f_smooth = traits.Either(float_trait, traits.List(float_trait),
                             default=2.0, argstr='-S %0.2f',
                             desc='smoothing kernel')
    f_decile = traits.Either(float_trait, traits.List(float_trait),
                             default=1.0, argstr='-d %0.5f',
                             desc=('set (decile) threshold to consider a '
                                   'computed gradient as outlier (ranges '
                                   '0.0-0.5)'))

    # Optimizer options
    iterations = traits.Either(traits.Int(), traits.List(traits.Int()),
                               default=50, argstr='-i %d',
                               desc='number of iterations (per level)')
    convergence_window = traits.Either(int_trait, traits.List(int_trait),
                                       default=10, argstr='-w %d',
                                       desc='convergence window in iterations')
    convergence_energy = traits.Either(bool_trait, traits.List(bool_trait),
                                       default=False,
                                       argstr='--convergence-energy',
                                       desc=('use lazy (default) or rigurous '
                                             'energy tracking'))
    grid_size = traits.Either(int_trait, traits.List(int_trait), default=8,
                              desc=('bspline control points per dimension and '
                                    'level'), argstr='-g %d')
    descript_update = traits.Either(int_trait, traits.List(int_trait),
                                    argstr='-u %d',
                                    desc=('update descriptors every N '
                                          'iterations, per level'))
    step_size = traits.Either(float_trait, traits.List(float_trait),
                              default=1.0, argstr='-s %0.5f',
                              desc=('update step size in gradient descent '
                                    'optimization'))
    alpha = traits.Either(float_trait, traits.List(float_trait), default=1.0,
                          argstr='-a %0.5f', desc='alpha scalar')
    beta = traits.Either(float_trait, traits.List(float_trait), default=1.0,
                         argstr='-b %0.5f', desc='beta scalar')


class ACWERegInputSpec(ACWERegInputGroupSpec):
    in_fixed = InputMultiPath(File(exists=True), argstr="-F %s",
                              mandatory=True,
                              desc=('target volume/image(s) contrast to '
                                    'register contours to'))
    in_prior = InputMultiPath(File(exists=True), argstr="-M %s",
                              mandatory=True,
                              desc=('vtk contours that will be registered to '
                                    'in_fixed. Should be given in hierarchical'
                                    ' order (from top to bottom, last is bg)'))
    levels = traits.Int(1, desc='number of levels in multi-resolution \
                        schemes', argstr="-L %d")
    out_prefix = traits.Str("regseg", desc='output files prefix',
                            argstr="-o %s", usedefault=True)
    log_filename = File(desc='filepath for log file', argstr='-l %s')
    images_verbosity = traits.Int(1, argstr='-v %d',
                                  desc=('verbosity of intermediate results '
                                        'output'))


class ACWERegOutputSpec(TraitedSpec):
    out_warped = OutputMultiPath(File(exists=True,
                                      desc='source images unwarped'))
    out_surfs = OutputMultiPath(File(exists=True,
                                     desc='priors in target space'))
    out_tpms = OutputMultiPath(File(exists=True,
                                    desc=('tissue probability maps (TPM) in '
                                          'target space')))
    out_field = File(exists=True, desc='output field')
    out_log = File(exists=True, desc='log JSON file')
    out_coeff = OutputMultiPath(File(desc='output coefficients'))


class ACWEReg(ANTSCommand):

    """
    Wraps regseg application from ACWERegistration to perform
    joint segmentation-registration based on ACWE.

    Example
    -------
    >>> regseg = ACWEReg()
    >>>
    'regseg -F T1w.nii.gz T2w.nii.gz -M csf.vtk white_lh.vtk white_rh.vtk \
pial_lh.vtk pial_rh.vtk -o tests [ -i 30 -u 10 -f 1.0 -s 0.5 -a 0.0 -b 0.0 \
-g 8 ]'
    """
    input_spec = ACWERegInputSpec
    input_group_spec = ACWERegInputGroupSpec
    output_spec = ACWERegOutputSpec
    _grouped_traits = []
    _cmd = 'regseg'
    _num_levels = 0

    def __init__(self, command=None, **inputs):
        """ Combine general and grouped inputs """
        super(ACWEReg, self).__init__(command=command, **inputs)
        self.groups = self.input_group_spec()

        general_names = ANTSCommandInputSpec().trait_names()

        for name in self.groups.trait_names():
            if not any(name in s for s in general_names):
                self._grouped_traits.append(name)

    def _parse_inputs(self, skip=None):
        """
        Parse all inputs using the ``argstr`` format string in the Trait.

        Any inputs that are assigned (not the default_value) are formatted
        to be added to the command line.

        Returns
        -------
        all_args : list
            A list of all inputs formatted for the command line.

        """
        all_args = []

        if skip is None:
            skip = []

        if isdefined(self.inputs.levels):
            self._num_levels = self.inputs.levels
        elif isdefined(self.inputs.iterations):
            if isinstance(self.inputs.iterations, list):
                self._num_levels = len(self.inputs.iterations)
            elif type(self.inputs.iterations) is int:
                self._num_levels = 1
            else:
                raise RuntimeError('iterations is not a valid value')
        else:
            raise RuntimeError('No way to guess number of levels')

        skip += ['levels']

        all_args += super(ACWEReg,
                          self)._parse_inputs(skip=skip + self._grouped_traits)

        for i in range(self._num_levels):
            all_args += [self._parse_group(i)]

        return all_args

    def _parse_group(self, gid, skip=None):
        retval = []
        retval.append(' [')

        metadata = dict(argstr=lambda t: t is not None)
        for name, spec in sorted(self.inputs.traits(**metadata).items()):
            if skip and name in skip:
                continue
            if name not in self._grouped_traits:
                continue

            value = getattr(self.inputs, name)

            if not isdefined(value):
                continue

            value = np.atleast_1d(value)

            if not len(value) == self._num_levels:
                raise RuntimeError(
                    'spec  \'%s\' should match number of levels' % name)

            val = value[gid]
            argval = self._format_group_arg(name, spec, val)
            retval.append(' ' + argval)

        retval.append(']')
        return "".join(retval)

    def _format_group_arg(self, name, spec, value):
        if isinstance(value, bool) or isinstance(value, np.bool_):
            if value:
                return spec.argstr
            else:
                return ''

        if value is None or value.__class__.__name__ == 'NoneType':
            return ''

        return super(ACWEReg, self)._format_arg(name, spec, value)

    def _format_arg(self, name, spec, value):
        if name == 'in_prior':
            idxs = []
            for i, v in enumerate(value):
                if 'white' in op.basename(v) or 'wm' in op.basename(v):
                    idxs.append(i)

            for i in reversed(idxs):
                value.append(value.pop(i))

            idxs = []
            for i, v in enumerate(value):
                if 'pial' in op.basename(v):
                    idxs.append(i)

            for i in reversed(idxs):
                value.append(value.pop(i))
        return super(ACWEReg, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        out_prefix = self.inputs.out_prefix
        outputs = self.output_spec().get()
        outputs['out_warped'] = [op.abspath('%s_warped_%d.nii.gz' % (
            out_prefix, i)) for i in range(len(self.inputs.in_fixed))]
        outputs['out_tpms'] = [op.abspath('%s_final_tpm_%d.nii.gz' % (
            out_prefix, i)) for i in range(len(self.inputs.in_prior) + 1)]
        outputs['out_surfs'] = [op.abspath('%s_warped_%s' % (
            out_prefix, op.basename(name))) for name in self.inputs.in_prior]
        outputs['out_field'] = op.abspath('%s_field.nii.gz' % out_prefix)
        outputs['out_coeff'] = op.abspath('%s_coeff.nii.gz' % out_prefix)

        logname = ''
        if isdefined(self.inputs.log_filename):
            logname = self.inputs.log_filename

        outputs['out_log'] = op.abspath('%s%s.log' % (out_prefix, logname))
        return outputs
