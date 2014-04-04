#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: Oscar Esteban - code@oscaresteban.es
# @Date:   2014-03-12 13:20:04
# @Last Modified by:   oesteban
# @Last Modified time: 2014-03-31 13:36:36

import os
import os.path as op
from glob import glob
import warnings

import numpy as np
import nibabel as nib

from nipype.interfaces.base import (traits, TraitedSpec, CommandLine,
                                    CommandLineInputSpec, InputMultiPath,
                                    OutputMultiPath, File,
                                    isdefined, Undefined )

from nipype.utils.filemanip import load_json, save_json, split_filename, fname_presuffix

warn = warnings.warn
warnings.filterwarnings('always', category=UserWarning)

class ACWERegInputGroupSpec( CommandLineInputSpec ):
    # Functional options
    float_trait = traits.Either( None, traits.Float(1.0) )
    int_trait = traits.Either( None, traits.Int( 0 ) )
    bool_trait = traits.Either( None, traits.Bool(False) )

    f_scale = traits.Either( float_trait, traits.List(float_trait), default=1.0,
                                desc='scales to be applied to computed shape gradients',
                                argstr='-f %0.5f')
    f_smooth = traits.Either( float_trait, traits.List(float_trait), default=2.0,
                                 desc='smoothing kernel', argstr='-S %0.2f' )
    f_decile = traits.Either( float_trait, traits.List(float_trait), default=1.0,
                                desc='set (decile) threshold to consider a computed \
                                      gradient as outlier (ranges 0.0-0.5)',
                                argstr='-d %0.5f')

    # Optimizer options
    iterations = traits.Either( traits.Int(), traits.List(traits.Int()), default=50,
                                desc='number of iterations (per level)', argstr='-i %d' )
    convergence_window = traits.Either( int_trait, traits.List( int_trait ), default=10,
                                desc='convergence window in iterations', argstr='-w %d' )
    convergence_energy = traits.Either( bool_trait, traits.List( bool_trait ), default=False,
                                        desc='use lazy (default) or rigurous energy tracking',
                                        argstr='--convergence-energy')
    grid_size = traits.Either( int_trait, traits.List( int_trait ), default=8,
                                desc='bspline control points per dimension and level',
                                argstr='-g %d')
    descript_update = traits.Either( int_trait, traits.List( int_trait ), default=5,
                                desc='update descriptors every N iterations, per level',
                                argstr='-u %d')
    step_size = traits.Either( float_trait, traits.List(float_trait), default=1.0,
                                desc='update step size in gradient descent optimization',
                                argstr='-s %0.5f')
    alpha = traits.Either( float_trait, traits.List(float_trait), default=1.0,
                                desc='alpha scalar',
                                argstr='-a %0.5f')
    beta = traits.Either( float_trait, traits.List(float_trait), default=1.0,
                                desc='beta scalar',
                                argstr='-b %0.5f')

class ACWERegInputSpec( ACWERegInputGroupSpec ):
    in_fixed = InputMultiPath(File(exists=True), argstr="-F %s", mandatory=True,
              desc='target volume/image(s) contrast to register contours to')
    in_prior = InputMultiPath(File(exists=True), argstr="-M %s", mandatory=True,
              desc='vtk contours that will be registered to in_fixed, should be \
                    given in hierarchical order (from top to bottom, last is bg)')
    levels = traits.Int(1, desc='number of levels in multi-resolution \
                        schemes', argstr="-L %d")
    out_prefix = traits.Str( "regseg", desc='output files prefix', argstr="-o %s",
         usedefault=True )
    log_filename = File(desc='filepath for log file', argstr='-l %s')
    images_verbosity = traits.Int(1, desc='verbosity of intermediate results output', argstr='-v %d')



class ACWERegOutputSpec( TraitedSpec ):
    out_warped = OutputMultiPath(File(exists=True, desc='source images unwarped'))
    out_surfs = OutputMultiPath(File(exists=True, desc='priors in target space'))
    out_tpms = OutputMultiPath(File(exists=True, desc='tissue probability maps (TPM) in target space'))
    out_field = File(exists=True, desc='output field')
    out_log = File(exists=True, desc='log JSON file')
    #out_coeff = OutputMultiPath( File(desc='output coefficients') )



class ACWEReg( CommandLine ):
    """ Wraps regseg application from ACWERegistration to perform
    joint segmentation-registration based on ACWE.

    Example
    -------
    >>> regseg = ACWEReg()
    >>>
    'regseg -F T1w.nii.gz T2w.nii.gz -M csf.vtk white_lh.vtk white_rh.vtk pial_lh.vtk pial_rh.vtk -o tests [ -i 30 -u 10 -f 1.0 -s 0.5 -a 0.0 -b 0.0 -g 8 ]'
    """
    input_spec = ACWERegInputSpec
    input_group_spec = ACWERegInputGroupSpec
    output_spec = ACWERegOutputSpec
    _grouped_traits = []
    _cmd = 'regseg'
    _num_levels = 0

    def __init__( self, command=None, **inputs ):
        """ Combine general and grouped inputs """
        super( ACWEReg, self ).__init__( command=command, **inputs )
        self.groups = self.input_group_spec()

        general_names = CommandLineInputSpec().trait_names()

        for name in self.groups.trait_names():
            if not any( name in s for s in general_names ):
                self._grouped_traits.append( name )


    def _parse_inputs(self, skip=None):
        """Parse all inputs using the ``argstr`` format string in the Trait.

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

        if isdefined( self.inputs.levels ):
            self._num_levels = self.inputs.levels
        elif isdefined( self.inputs.iterations ):
            if isinstance( self.inputs.iterations, list ):
                self._num_levels = len( self.inputs.iterations )
            elif type( self.inputs.iterations ) is int:
                self._num_levels = 1
            else:
                raise RuntimeError( 'iterations is not a valid value')
        else:
            raise RuntimeError( 'No way to guess number of levels')

        skip+=['levels']

        all_args+=super( ACWEReg, self )._parse_inputs( skip=skip+self._grouped_traits )

        for i in range( self._num_levels ):
            all_args+=[ self._parse_group( i ) ]

        return all_args


    def _parse_group( self, gid, skip=None ):
        retval = []
        retval.append( ' [' )

        metadata = dict(argstr=lambda t: t is not None)
        for name, spec in sorted(self.inputs.traits(**metadata).items()):
            if skip and name in skip:
                continue
            if not name in self._grouped_traits:
                continue

            value = getattr( self.inputs, name)

            if not isdefined(value):
                continue

            value = np.atleast_1d( value )

            if not len(value) == self._num_levels:
                raise RuntimeError('spec  \'%s\' should match number of levels' % name )

            val = value[gid]
            argval = self._format_group_arg( name, spec, val )
            retval.append( ' ' +  argval )

        retval.append( ']' )
        return "".join(retval)

    def _format_group_arg( self, name, spec, value ):
        if isinstance(value,bool) or isinstance(value,np.bool_):
            if value:
                return spec.argstr
            else:
                return ''

        if value is None or value.__class__.__name__=='NoneType':
            return ''

        return super( ACWEReg, self )._format_arg( name, spec, value )

    def _format_arg( self, name, spec, value ):
        if name == 'in_prior':
            idxs = []
            for i,v in enumerate(value):
                if 'white' in op.basename(v) or 'wm' in op.basename(v):
                    idxs.append( i )

            for i in reversed(idxs):
                value.append( value.pop( i ) )

            idxs = []
            for i,v in enumerate(value):
                if 'pial' in op.basename(v):
                    idxs.append( i )

            for i in reversed(idxs):
                value.append( value.pop( i ) )
        return super( ACWEReg, self )._format_arg( name, spec, value )

    def _list_outputs(self):
        out_prefix = self.inputs.out_prefix
        outputs = self.output_spec().get()
        outputs['out_warped'] = [ op.abspath('%s_warped_%d.nii.gz' % (out_prefix, i))  for i in range(len(self.inputs.in_fixed )) ]
        outputs['out_tpms'] = [ op.abspath('%s_final_tpm_%d.nii.gz' % (out_prefix, i))  for i in range(len(self.inputs.in_prior) + 1 ) ]
        outputs['out_surfs'] = [ op.abspath('%s_warped_%s' % (out_prefix, op.basename(name)))  for name in self.inputs.in_prior ]
        outputs['out_field'] = op.abspath('%s_displacement_field.nii.gz' % out_prefix )

        logname = ''
        if isdefined( self.inputs.log_filename ):
            logname = self.inputs.log_filename

        outputs['out_log'] = op.abspath( '%s%s.log' % (out_prefix, logname) )

        return outputs


class RandomBSplineDeformationInputSpec( CommandLineInputSpec ):
    in_file = InputMultiPath(File(exists=True), mandatory=True, argstr="-I %s",
              desc='image(s) to be deformed')
    in_coeff = InputMultiPath(File(exists=True), mandatory=False, argstr="-C %s",
               desc='coefficient images (to re-use a deformation field)')
    in_surfs = InputMultiPath(File(exists=True), mandatory=False, argstr="-S %s",
               desc='apply deformation field to surfaces')
    in_mask = File(exists=True, argstr='-M %s', desc='set a mask')
    grid_size_item_trait = traits.Int(10, usedefault=True )
    grid_size = traits.Either( grid_size_item_trait, traits.List(grid_size_item_trait),
                               xor=['in_coeff'], default=10, usedefault=True,
                               desc='size of control points grid', argstr="-g %s" )
    out_prefix = traits.Str( "def", desc='output files prefix', argstr="-o %s", usedefault=True,
                             mandatory=True )

class RandomBSplineDeformationOutputSpec( TraitedSpec ):
    out_file = OutputMultiPath(File(exists=True, desc='warped input files'))
    out_coeff = OutputMultiPath(File(exists=True, desc='output coefficients'))
    out_field = File(exists=True, desc='output warping field')
    out_field_base = OutputMultiPath(File(exists=True, desc='output warping bspline field'))
    out_surfs = OutputMultiPath(File(desc='output warped surfaces'))
    out_mask = File(exists=True, desc='warped input mask')

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
    _cmd = 'bspline_field'

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
        outputs['out_field'] = op.abspath( '%s_dispfield.nii.gz' % out_prefix )
        outputs['out_coeff'] = [ op.abspath( '%s_coeffs_%d.nii.gz' % ( out_prefix, i ))  for i in range(3) ]
        outputs['out_field_base'] = [ op.abspath( '%s_field_cmp%d.nii.gz' % ( out_prefix, i ))  for i in range(3) ]

        if isdefined( self.inputs.in_surfs ):
            outputs['out_surfs'] = [ op.abspath( '%s_surf_%d.vtk' % ( out_prefix, i ))  for i in range(len(self.inputs.in_surfs)) ]

        if isdefined( self.inputs.in_mask ):
            outputs['out_mask'] = op.abspath( '%s_mask_warped.nii.gz' % out_prefix )

        return outputs


class FieldBasedWarpInputSpec( CommandLineInputSpec ):
    in_file = InputMultiPath(File(exists=True), mandatory=True, argstr="-I %s",
              desc='image(s) to be deformed')
    in_field = File(exists=True, mandatory=False, argstr="-F %s",
               desc='field')
    in_mask = File(exists=True, argstr='-M %s', desc='set a mask')
    in_surf = InputMultiPath(File(exists=True), argstr="-S %s",
              desc='surface(s) to be deformed')
    out_prefix = traits.Str( "fbased", desc='output files prefix', argstr="-o %s", usedefault=True,
                             mandatory=True )

class FieldBasedWarpOutputSpec( TraitedSpec ):
    out_file = OutputMultiPath(File(exists=True, desc='warped input files'))
    out_surf = OutputMultiPath(File(exists=True, desc='warped input surfaces'))
    out_mask = File(exists=True, desc='warped input mask')

class FieldBasedWarp( CommandLine ):
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
    'warp_image -I moving.nii -f field.nii -S lh.white.vtk rh.white.vtk -o myprefix'
    """

    input_spec = FieldBasedWarpInputSpec
    output_spec = FieldBasedWarpOutputSpec
    _cmd = 'warp_image'

    def _list_outputs( self ):
        out_prefix = self.inputs.out_prefix
        outputs = self.output_spec().get()

        outputs['out_file'] = [ op.abspath( '%s_warped_%d.nii.gz' % ( out_prefix, i ))  for i in range(len(self.inputs.in_file)) ]

        if isdefined( self.inputs.in_mask ):
            outputs['out_mask'] = op.abspath( '%s_mask_warped.nii.gz' % out_prefix )

        if isdefined( self.inputs.in_surf ):
            outputs['out_surf'] = [ op.abspath( '%s_warped_%d.vtk' % ( out_prefix, i ))  for i in range(len(self.inputs.in_surf)) ]
        return outputs

