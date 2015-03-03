#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2015-01-15 10:47:12
# @Last Modified by:   oesteban
# @Last Modified time: 2015-03-03 14:37:13

import os.path as op


def hcp_workflow(name='Evaluation_HCP', settings={},
                 map_metric=False, compute_fmb=False):
    """
    The pyacwereg evaluation workflow for the human connectome project (HCP)
    """
    from nipype.pipeline import engine as pe
    from nipype.interfaces import utility as niu
    from nipype.interfaces import io as nio
    from nipype.interfaces import freesurfer as fs
    from nipype.algorithms.mesh import ComputeMeshWarp, WarpPoints
    from nipype.algorithms.misc import AddCSVRow
    from nipype.workflows.dmri.fsl.artifacts import sdc_fmb

    from pyacwereg import misc as acwregmisc
    from pyacwereg import data
    from pyacwereg.interfaces.utility import ExportSlices
    from pyacwereg.workflows.registration import regseg_wf, sdc_t2b
    from pyacwereg.workflows import evaluation as ev
    from pyacwereg.workflows.preprocess import preprocess

    from pysdcev.workflows.warpings import process_vsm
    from pysdcev.workflows.tractography import mrtrix_dti

    wf = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['subject_id', 'data_dir']), name='inputnode')
    inputnode.inputs.data_dir = settings['data_dir']
    inputnode.iterables = [('subject_id', settings['subject_id'])]

    # Generate the distorted set, including surfaces
    pre = preprocess()
    rdti = mrtrix_dti('ReferenceDTI')
    wdti = mrtrix_dti('WarpedDTI')
    mdti = pe.Node(niu.Merge(2), name='MergeDTI')

    wf.connect([
        (inputnode, pre, [('subject_id', 'inputnode.subject_id'),
                          ('data_dir', 'inputnode.data_dir')]),
        (pre,      rdti, [('outputnode.dwi', 'inputnode.in_dwi'),
                          ('outputnode.dwi_mask', 'inputnode.in_mask'),
                          ('outputnode.bvec', 'inputnode.in_bvec'),
                          ('outputnode.bval', 'inputnode.in_bval')]),
        (pre,      wdti, [('outputnode.warped_dwi', 'inputnode.in_dwi'),
                          ('outputnode.warped_msk', 'inputnode.in_mask'),
                          ('outputnode.bvec', 'inputnode.in_bvec'),
                          ('outputnode.bval', 'inputnode.in_bval')]),
        (wdti,     mdti, [('outputnode.fa', 'in1'),
                          ('outputnode.md', 'in2')]),
    ])

    regseg = regseg_wf(usemask=True)
    regseg.inputs.inputnode.options = data.get('regseg_hcp')
    exprs = pe.Node(ExportSlices(all_axis=True), name='ExportREGSEG')
    meshrs = pe.MapNode(ComputeMeshWarp(),
                        iterfield=['surface1', 'surface2'],
                        name='REGSEGSurfDistance')
    csvrs = pe.Node(AddCSVRow(in_file=settings['out_csv']),
                    name="REGSEGAddRow")
    csvrs.inputs.method = 'REGSEG'

    wf.connect([
        (mdti,      regseg, [('out', 'inputnode.in_fixed')]),
        (pre,       regseg, [('outputnode.surf', 'inputnode.in_surf'),
                             ('outputnode.warped_msk', 'inputnode.in_mask')]),
        (pre,        exprs, [('outputnode.warped_surf', 'surfaces0')]),
        (regseg,     exprs, [('outputnode.out_surf', 'surfaces1')]),
        (wdti,       exprs, [('outputnode.fa', 'reference')]),
        (pre,       meshrs, [('outputnode.warped_surf', 'surface1')]),
        (regseg,    meshrs, [('outputnode.out_surf', 'surface2')]),
        (inputnode,  csvrs, [('subject_id', 'subject_id')]),
        (meshrs,     csvrs, [('distance', 'surf_dist')])
    ])

    if compute_fmb:
        cmethod0 = sdc_fmb()
        selbmap = pe.Node(niu.Split(splits=[1, 1], squeeze=True),
                          name='SelectBmap')
        dfm = process_vsm()
        dfm.inputs.inputnode.scaling = 1.0
        dfm.inputs.inputnode.enc_dir = 'y-'
        wrpsurf = pe.MapNode(WarpPoints(), iterfield=['points'],
                             name='UnwarpSurfs')
        export0 = pe.Node(ExportSlices(all_axis=True), name='ExportFMB')
        mesh0 = pe.MapNode(ComputeMeshWarp(),
                           iterfield=['surface1', 'surface2'],
                           name='FMBSurfDistance')
        csv0 = pe.Node(AddCSVRow(in_file=settings['out_csv']),
                       name="FMBAddRow")
        csv0.inputs.method = 'FMB'

        wf.connect([
            (pre,       cmethod0, [
                ('outputnode.warped_dwi', 'inputnode.in_file'),
                ('outputnode.warped_msk', 'inputnode.in_mask'),
                ('outputnode.bval', 'inputnode.in_bval'),
                ('outputnode.mr_param', 'inputnode.settings')]),
            (pre,        selbmap, [('outputnode.bmap_wrapped', 'inlist')]),
            (selbmap,   cmethod0, [('out1', 'inputnode.bmap_mag'),
                                   ('out2', 'inputnode.bmap_pha')]),
            (cmethod0,       dfm, [('outputnode.out_vsm', 'inputnode.vsm')]),
            (pre,            dfm, [
                ('outputnode.warped_msk', 'inputnode.reference')]),
            (dfm,        wrpsurf, [('outputnode.dfm', 'warp')]),
            (pre,        wrpsurf, [('outputnode.surf', 'points')])
            (wrpsurf,    export0, [('out_points', 'surfaces0')]),
            (pre,        export0, [('outputnode.warped_surf', 'surfaces1')]),
            (wdti,       export0, [('outputnode.fa', 'reference')]),
            (pre,          mesh0, [('outputnode.warped_surf', 'surface1')]),
            (wrpsurf,      mesh0, [('out_points', 'surface2')]),
            (inputnode,     csv0, [('subject_id', 'subject_id')]),
            (mesh0,         csv0, [('distance', 'surf_dist')])
        ])

    cmethod1 = sdc_t2b(num_threads=settings['nthreads'])
    export1 = pe.Node(ExportSlices(all_axis=True), name='ExportT2B')
    mesh1 = pe.MapNode(ComputeMeshWarp(),
                       iterfield=['surface1', 'surface2'],
                       name='T2BSurfDistance')
    csv1 = pe.Node(AddCSVRow(in_file=settings['out_csv']),
                   name="T2BAddRow")
    csv1.inputs.method = 'T2B'

    wf.connect([
        (pre,       cmethod1, [
            ('outputnode.warped_dwi', 'inputnode.in_dwi'),
            ('outputnode.warped_msk', 'inputnode.dwi_mask'),
            ('outputnode.t2w_brain', 'inputnode.in_t2w'),
            ('outputnode.t1w_mask', 'inputnode.t2w_mask'),
            ('outputnode.surf', 'inputnode.in_surf'),
            ('outputnode.bval', 'inputnode.in_bval'),
            ('outputnode.mr_param', 'inputnode.in_param')]),
        (cmethod1,   export1, [('outputnode.out_surf', 'surfaces0')]),
        (pre,        export1, [('outputnode.warped_surf', 'surfaces1')]),
        (wdti,       export1, [('outputnode.fa', 'reference')]),
        (pre,          mesh1, [('outputnode.warped_surf', 'surface1')]),
        (cmethod1,     mesh1, [('outputnode.out_surf', 'surface2')]),
        (inputnode,     csv1, [('subject_id', 'subject_id')]),
        (mesh1,         csv1, [('distance', 'surf_dist')])
    ])

    if map_metric:
        out_csv = op.abspath(op.join(name, 'energiesmapping.csv'))
        mapen = ev.map_energy(out_csv=out_csv)
        wf.connect([
            (inputnode, mapen, [('subject_id', 'inputnode.subject_id')]),
            (regseg,    mapen, [('outputnode.out_enh', 'inputnode.reference'),
                                ('outputnode.reg_msk', 'inputnode.in_mask')]),
            (pre,       mapen, [
                ('outputnode.warped_surf', 'inputnode.surfaces0'),
                ('outputnode.surf', 'inputnode.surfaces1')])
        ])

    return wf
