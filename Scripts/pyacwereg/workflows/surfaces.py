#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: Oscar Esteban - code@oscaresteban.es
# @Date:   2014-03-05 15:08:55
# @Last Modified by:   Oscar Esteban
# @Last Modified time: 2015-03-03 15:40:11
"""
Surface extraction
++++++++++++++++++

Defines the workflows for extracting surfaces from segmentations

:platform: Unix
:moduleauthor: Oscar Esteban <code@oscaresteban>

"""

import os
import os.path as op
import nipype.pipeline.engine as pe             # pipeline engine
from nipype.interfaces import utility as niu    # utility
from nipype.interfaces import fsl               # fsl
from nipype.interfaces import freesurfer as fs  # freesurfer
from nipype.interfaces import ants              # ANTS


def extract_surface(name='GenSurface'):
    """ A nipype workflow for surface extraction from ``labels`` in a segmentation.

    .. note :: References used to implement this code:

      * <https://github.com/nipy/nipype/issues/307>
      * <https://mail.nmr.mgh.harvard.edu/pipermail/\
freesurfer/2011-November/021391.html>
      * <http://brainder.org/2012/05/08/importing-\
freesurfer-subcortical-structures-into-blender/>
      * <https://mail.nmr.mgh.harvard.edu/pipermail/\
freesurfer/2013-June/030586.html>
    """
    import nipype.pipeline.engine as pe
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['aseg', 'norm', 'in_filled', 'labels', 'name'],
        mandatory_inputs=False),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_surf', 'out_binary']), name='outputnode')
    binarize = pe.Node(fs.Binarize(), name='BinarizeLabels')
    fill = pe.Node(niu.Function(
        function=_fillmask, input_names=['in_file', 'in_filled'],
        output_names=['out_file']), name='FillMask')
    pretess = pe.Node(fs.MRIPretess(label=1), name='PreTess')
    tess = pe.Node(fs.MRITessellate(label_value=1), name='tess')
    smooth = pe.Node(fs.SmoothTessellation(disable_estimates=True),
                     name='mris_smooth')
    rename = pe.Node(niu.Rename(keep_ext=False),
                     name='rename')

    tovtk = pe.Node(fs.MRIsConvert(out_datatype='vtk'), name='toVTK')
    fixVTK = pe.Node(niu.Function(
        input_names=['in_file', 'in_ref'], output_names=['out_file'],
        function=_fixvtk), name='fixVTK')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,   binarize,   [('aseg', 'in_file'),
                                   ('labels', 'match')]),
        (inputnode,   fixVTK,     [('norm', 'in_ref')]),
        (inputnode,   pretess,    [('norm', 'in_norm')]),
        (inputnode,   fill,       [('in_filled', 'in_filled')]),
        (binarize,    fill,       [('binary_file', 'in_file')]),
        (fill,        pretess,    [('out_file', 'in_filled')]),
        (pretess,     tess,       [('out_file', 'in_file')]),
        (tess,        smooth,     [('surface', 'in_file')]),
        (smooth,      rename,     [('surface', 'in_file')]),
        (inputnode,   rename,     [('name', 'format_string')]),
        (rename,      tovtk,      [('out_file', 'in_file')]),
        (tovtk,       fixVTK,     [('converted', 'in_file')]),
        (fixVTK,      outputnode, [('out_file', 'out_surf')]),
        (fill,        outputnode, [('out_file', 'out_binary')])
    ])
    return wf


def all_surfaces(name='Surfaces', gen_outer=False):
    import nipype.pipeline.engine as pe
    from nipype.interfaces.io import JSONFileGrabber
    import pyacwereg.data as data

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['aseg', 'norm', 'in_mask']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_surf']), name='outputnode')

    readls = pe.Node(JSONFileGrabber(in_file=data.get('model_labels')),
                     name='ReadModelLabels')

    nsurfs = 0
    tha = extract_surface(name='ThalSurface')
    tha.inputs.inputnode.name = '%02d.thalamus' % nsurfs
    nsurfs += 1

    csf = extract_surface(name='DarkFASurface')
    csf.inputs.inputnode.name = '%02d.csf_dgm' % nsurfs
    nsurfs += 1

    bstem = extract_surface(name='stemSurface')
    bstem.inputs.inputnode.name = '%02d.bstem' % nsurfs
    nsurfs += 1

    wm = extract_surface(name='WMSurface')
    wm.inputs.inputnode.name = '%02d.white' % nsurfs
    nsurfs += 1

    cgm = extract_surface(name='cGMSurface')
    cgm.inputs.inputnode.name = '%02d.cgm' % nsurfs
    nsurfs += 1

    pial = extract_surface(name='PialSurface')
    pial.inputs.inputnode.name = '%02d.pial' % nsurfs
    nsurfs += 1

    if gen_outer:
        nsurfs = nsurfs + 1

    m = pe.Node(niu.Merge(nsurfs), name='MergeSurfs')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, tha,   [('aseg', 'inputnode.aseg'),
                            ('norm', 'inputnode.norm')]),
        (readls,    tha,   [('thal_labels', 'inputnode.labels')]),
        (tha,       m,     [('outputnode.out_surf', 'in1')]),
        (inputnode, csf,   [('aseg', 'inputnode.aseg'),
                            ('norm', 'inputnode.norm')]),
        (readls,    csf,   [('csf_dgm_labels', 'inputnode.labels')]),
        (csf,       m,     [('outputnode.out_surf', 'in2')]),
        (inputnode, bstem, [('aseg', 'inputnode.aseg'),
                            ('norm', 'inputnode.norm')]),
        (readls,    bstem, [('bstem_labels', 'inputnode.labels')]),
        (bstem,     m,     [('outputnode.out_surf', 'in3')]),
        (inputnode, wm,    [('aseg', 'inputnode.aseg'),
                            ('norm', 'inputnode.norm')]),
        (readls,    wm,    [('wm_labels', 'inputnode.labels')]),
        (wm,        m,     [('outputnode.out_surf', 'in4')]),
        (inputnode, cgm,   [('aseg', 'inputnode.aseg'),
                            ('norm', 'inputnode.norm')]),
        (readls,    cgm,   [('cgm_labels', 'inputnode.labels')]),
        (cgm,       m,     [('outputnode.out_surf', 'in5')]),
        (inputnode, pial,  [('aseg', 'inputnode.aseg'),
                            ('norm', 'inputnode.norm')]),
        (readls,    pial,  [('gm_labels', 'inputnode.labels')]),
        (pial,      m,     [('outputnode.out_surf', 'in6')]),
        (m,    outputnode, [('out', 'out_surf')])
    ])

    if gen_outer:
        msk = extract_surface(name='MaskSurf')
        msk.inputs.inputnode.labels = [1]
        msk.inputs.inputnode.name = '%01d.outer' % nsurfs - 1

        wf.connect([
            (inputnode, msk,  [('in_mask', 'inputnode.aseg'),
                               ('in_mask', 'inputnode.norm')]),
            (msk,       m,    [('outputnode.out_surf', 'in%d' % nsurfs)])
        ])
    return wf


def _fillmask(in_file, in_filled=None):
    import nibabel as nb
    import numpy as np
    from nipype.interfaces.base import isdefined
    import os.path as op

    if (not isdefined(in_filled)) or (in_filled is None):
        return in_file

    nii = nb.load(in_file)
    data = nii.get_data()

    in_filled = np.atleast_1d(in_filled).tolist()
    for fname in in_filled:
        data = data + nb.load(fname).get_data()
    data[data > 1.0] = 1.0

    out_file = op.abspath('mask_filled.nii.gz')
    nb.Nifti1Image(data.astype(np.uint8), nii.get_affine(),
                   nii.get_header()).to_filename(out_file)
    return out_file


def _fixvtk(in_file, in_ref, out_file=None):
    """
    Transforms a vtk file from Freesurfer's *tkRAS* coordinates
    to a target image in *scannerRAS* coordinates.
    """
    import nibabel as nb
    import numpy as np
    import os.path as op
    import subprocess as sp

    if out_file is None:
        fname, ext = op.splitext(op.basename(in_file))
        if ext == ".gz":
            fname, _ = op.splitext(fname)

        out_file = op.abspath("%s_fixed.vtk" % fname)

    ref = nb.load(in_ref)
    cmd_info = "mri_info %s  --tkr2scanner" % in_ref
    proc = sp.Popen(cmd_info, stdout=sp.PIPE, shell=True)
    data = bytearray(proc.stdout.read())
    if 'niiRead' in data:
        _, data = data.split('\n', 1)

    mstring = np.fromstring(data.decode("utf-8"), sep='\n')
    matrix = np.reshape(mstring, (4, -1))

    with open(in_file, 'r') as f:
        with open(out_file, 'w+') as w:
            npoints = 0
            pointid = -5

            for i, l in enumerate(f):
                if (i == 4):
                    s = l.split()
                    npoints = int(s[1])
                    fmt = np.dtype(s[2])
                elif((i > 4) and (pointid < npoints)):
                    vert = [float(x) for x in l.split()]
                    vert.append(1.0)
                    newvert = np.dot(matrix, vert)
                    l = '%.9f  %.9f  %.9f\n' % tuple(newvert[0:3])

                w.write(l)
                pointid = pointid + 1

    return out_file
