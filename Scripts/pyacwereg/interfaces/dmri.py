#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2015-03-10 16:15:07
# @Last Modified by:   oesteban
# @Last Modified time: 2015-03-11 10:30:33

import os
import os.path as op
import nibabel as nb
import numpy as np

from nipype.interfaces.base import (BaseInterface, traits, TraitedSpec, File,
                                    InputMultiPath, OutputMultiPath,
                                    BaseInterfaceInputSpec, isdefined,
                                    DynamicTraitedSpec, Directory,
                                    CommandLine, CommandLineInputSpec)

from nipype import logging
iflogger = logging.getLogger('interface')


class PhaseUnwrapInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True,
                   desc='phase file to be unwrapped')
    in_mask = File(exists=True, desc='mask file')
    rescale = traits.Bool(True, usedefault=True,
                          desc='rescale range to 2*pi')
    out_file = File('unwrapped.nii.gz', usedefault=True,
                    desc='output file name')


class PhaseUnwrapOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='output file name')


class PhaseUnwrap(BaseInterface):

    """
    Unwraps a phase file
    """
    input_spec = PhaseUnwrapInputSpec
    output_spec = PhaseUnwrapOutputSpec

    def _run_interface(self, runtime):
        from skimage.restoration import unwrap_phase as unwrap
        from skimage.restoration import nl_means_denoising as denoise
        from math import pi

        im = nb.load(self.inputs.in_file)
        wrapped = im.get_data()

        if self.inputs.rescale:
            wrapped = wrapped - wrapped.min()
            wrapped = (2.0 * pi * wrapped) / wrapped.max()
            wrapped -= pi

        unw = unwrap(wrapped).astype(np.float32)

        msk = None
        if isdefined(self.inputs.in_mask):
            msk = nb.load(self.inputs.in_mask).get_data()
            msk[msk > 0.0] = 1.0
            msk[msk < 1.0] = 0.0
            unw *= msk
            unw = np.ma.array(wrapped, mask=1-msk)

        unw = denoise(unw, 7, h=0.15, multichannel=False).astype(np.float32)

        # if msk is not None:
        #     unw = np.ma.array(unw, mask=np.zeros_like(msk))
        #     unw[msk < 1.0] = 0

        hdr = im.get_header().copy()
        hdr.set_data_dtype(np.float32)
        nb.Nifti1Image(unw.astype(np.float32), im.get_affine(),
                       hdr).to_filename(op.abspath(self.inputs.out_file))
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = op.abspath(self.inputs.out_file)
        return outputs
