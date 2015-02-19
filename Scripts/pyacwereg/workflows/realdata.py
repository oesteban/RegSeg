#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2015-01-15 10:47:12
# @Last Modified by:   Oscar Esteban
# @Last Modified time: 2015-02-19 19:17:37


def hcp_workflow(name='Evaluation_HCP', settings={}):
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
    from pyacwereg.interfaces.utility import ExportSlices
    from pyacwereg.workflows.registration import regseg_wf, sdc_t2b
    from pyacwereg.workflows import evaluation as ev
    from pyacwereg.workflows.fieldmap import bmap_registration

    from pysdcev.workflows.warpings import process_vsm
    from pysdcev.workflows.smri import preprocess_t2, preprocess_dwi
    from pysdcev.workflows.tractography import mrtrix_dti
    from pysdcev.stages.stage1 import stage1

    fnames = dict(t1w='T1w_acpc_dc_restore.nii.gz',
                  t1w_brain='T1w_acpc_dc_restore_brain.nii.gz',
                  t2w='T2w_acpc_dc_restore.nii.gz',
                  t2w_brain='T2w_acpc_dc_restore_brain.nii.gz',
                  t1w_mask='brainmask_fs.nii.gz',
                  dwi='dwi.nii.gz',
                  dwi_mask='dwi_brainmask.nii.gz',
                  bval='bvals',
                  bvec='bvecs',
                  aseg='aparc+aseg.nii.gz',
                  fmap_mag='FM_mag.nii.gz',
                  fmap_pha='FM_pha.nii.gz',
                  mr_param='parameters.txt')

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['subject_id', 'data_dir']), name='inputnode')
    inputnode.inputs.data_dir = settings['data_dir']
    inputnode.iterables = [('subject_id', settings['subject_id'])]

    ds_tpl_args = {k: [['subject_id', [v]]] for k, v in fnames.iteritems()}
    ds = pe.Node(nio.DataGrabber(
        infields=['subject_id'], outfields=ds_tpl_args.keys(),
        sort_filelist=False, template='*'), name='sMRISource')
    ds.inputs.field_template = {k: 'subjects/%s/%s'
                                for k in ds_tpl_args.keys()}
    ds.inputs.template_args = ds_tpl_args

    bmap_prep = bmap_registration()
    bmap_prep.inputs.inputnode.factor = 4.0

    poly_msk = pe.Node(fs.Binarize(), name='GenPolyMask')
    poly_msk.inputs.match = [4, 5, 43, 44, 14, 72, 24, 2, 28, 31,
                             41, 60, 63, 77, 78, 79, 85, 86, 100, 108, 109,
                             117, 250, 251, 252, 253, 254, 255,
                             9, 10, 48, 49, 107, 116,        # thalamus
                             13, 52, 104, 113,               # pallidum
                             12, 51, 102, 111, 136, 137,     # putamen
                             11, 50, 101, 110,               # caudate
                             26, 58, 103, 112,               # accumbens
                             17, 53, 106, 115,               # hippocampus
                             18, 54, 96, 97, 105, 114, 218,
                             30, 62, 72, 75, 76, 98,
                             ] + range(1000, 1036) + range(2000, 2036)
    pmregrid = pe.Node(fs.MRIConvert(out_type='niigz', out_datatype='uchar'),
                       name='RegridPolymask')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,  ds,        [('subject_id', 'subject_id'),
                                 ('data_dir', 'base_directory')]),
        (ds,         poly_msk,  [('aseg', 'in_file')]),
        (poly_msk,   pmregrid,  [('binary_file', 'in_file')]),
        (ds,         pmregrid,  [('dwi_mask', 'reslice_like')]),
        (ds,         bmap_prep, [('t1w_brain', 'inputnode.t1w_brain'),
                                 ('dwi_mask', 'inputnode.dwi_mask'),
                                 ('fmap_mag', 'inputnode.mag'),
                                 ('fmap_pha', 'inputnode.pha')]),
        (pmregrid,   bmap_prep, [('out_file', 'inputnode.poly_mask')])
    ])

    rdti = mrtrix_dti('ReferenceDTI')
    wf.connect([
        (ds,    rdti,   [('dwi', 'inputnode.in_dwi'),
                         ('dwi_mask', 'inputnode.in_mask'),
                         ('bvec', 'inputnode.in_bvec'),
                         ('bval', 'inputnode.in_bval')])
    ])

    st1 = stage1()
    wf.connect([
        (ds,        st1, [('t1w', 'inputnode.t1w'),
                          ('t2w', 'inputnode.t2w'),
                          ('t1w_brain', 'inputnode.t1w_brain'),
                          ('t2w_brain', 'inputnode.t2w_brain'),
                          ('t1w_mask', 'inputnode.t1w_mask'),
                          ('dwi', 'inputnode.dwi'),
                          ('dwi', 'inputnode.dwi_brain'),
                          ('dwi_mask', 'inputnode.dwi_mask'),
                          ('bvec', 'inputnode.bvec'),
                          ('bval', 'inputnode.bval'),
                          ('aseg', 'inputnode.aseg'),
                          ('aseg', 'inputnode.parcellation'),
                          ('mr_param', 'inputnode.mr_params')]),
        (bmap_prep, st1, [
            ('outputnode.wrapped', 'inputnode.bmap_wrapped'),
            ('outputnode.unwrapped', 'inputnode.bmap_unwrapped')])
    ])

    dti = mrtrix_dti()
    mdti = pe.Node(niu.Merge(2), name='MergeDTI')

    regseg = regseg_wf(usemask=True)
    regseg.inputs.inputnode.images_verbosity = 3
    regseg.inputs.inputnode.alpha = [0.0, 0.0, 0.0]
    regseg.inputs.inputnode.beta = [0.0, 0.0, 0.0]
    regseg.inputs.inputnode.convergence_energy = [True, True, True]
    regseg.inputs.inputnode.descript_adaptative = [True, True, False]
    regseg.inputs.inputnode.convergence_value = [1.e-6, 1.e-7, 1.e-8]
    regseg.inputs.inputnode.convergence_window = [50, 20, 5]
    regseg.inputs.inputnode.descript_update = [None, None, None]
    regseg.inputs.inputnode.f_smooth = [2.0, 0.5, None]
    regseg.inputs.inputnode.grid_spacing = [
        (40., 100., 40.), (30., 30., 30.), (20., 30., 10.)]
    regseg.inputs.inputnode.iterations = [350, 250, 100]
    regseg.inputs.inputnode.scales = [(0.0, 1.0, 0.0)] * 3
    regseg.inputs.inputnode.step_size = [1.e-5, 5.e-5, 1.e-4]

    wf.connect([
        (st1,    dti,    [('out_dis_set.dwi', 'inputnode.in_dwi'),
                          ('out_dis_set.dwi_mask', 'inputnode.in_mask')]),
        (ds,     dti,    [('bvec', 'inputnode.in_bvec'),
                          ('bval', 'inputnode.in_bval')]),
        (dti,    mdti,   [('outputnode.fa', 'in1'),
                          ('outputnode.md', 'in2')]),
        (mdti,   regseg, [('out', 'inputnode.in_fixed')]),
        (st1,    regseg, [('out_dis_set.tpms', 'inputnode.in_tpms'),
                          ('out_ref_set.surf', 'inputnode.in_surf'),
                          ('out_dis_set.dwi_mask', 'inputnode.in_mask')])
    ])

    cmethod0 = sdc_fmb()
    selbmap = pe.Node(niu.Split(splits=[1, 1], squeeze=True),
                      name='SelectBmap')
    dfm = process_vsm()
    dfm.inputs.inputnode.scaling = 1.0
    dfm.inputs.inputnode.enc_dir = 'y-'
    wrpsurf = pe.MapNode(WarpPoints(), iterfield=['points'],
                         name='UnwarpSurfs')

    wf.connect([
        (st1,       cmethod0, [('out_dis_set.dwi', 'inputnode.in_file'),
                               ('out_dis_set.dwi_mask', 'inputnode.in_mask')]),
        (ds,        cmethod0, [('bval', 'inputnode.in_bval'),
                               ('mr_param', 'inputnode.settings')]),
        (bmap_prep,  selbmap, [('outputnode.wrapped', 'inlist')]),
        (selbmap,   cmethod0, [('out1', 'inputnode.bmap_mag'),
                               ('out2', 'inputnode.bmap_pha')]),
        (cmethod0,       dfm, [('outputnode.out_vsm', 'inputnode.vsm')]),
        (st1,            dfm, [
            ('out_dis_set.dwi_mask', 'inputnode.reference')]),
        (dfm,        wrpsurf, [('outputnode.dfm', 'warp')]),
        (st1,        wrpsurf, [('out_ref_set.surf', 'points')])
    ])

    export0 = pe.Node(ExportSlices(all_axis=True), name='ExportREGSEG')
    export1 = pe.Node(ExportSlices(all_axis=True), name='ExportFMB')

    wf.connect([
        (regseg,   export0, [('outputnode.out_surf', 'surfaces0')]),
        (st1,      export0, [('out_dis_set.surf', 'surfaces1')]),
        (dti,      export0, [('outputnode.fa', 'reference')]),
        (wrpsurf,  export1, [('out_points', 'surfaces0')]),
        (st1,      export1, [('out_dis_set.surf', 'surfaces1')]),
        (dti,      export1, [('outputnode.fa', 'reference')])
    ])

    mesh0 = pe.MapNode(ComputeMeshWarp(),
                       iterfield=['surface1', 'surface2'],
                       name='REGSEGSurfDistance')
    csv0 = pe.Node(AddCSVRow(in_file=settings['out_csv']),
                   name="REGSEGAddRow")
    csv0.inputs.method = 'REGSEG'

    wf.connect([
        (st1,       mesh0, [('out_dis_set.surf', 'surface1')]),
        (regseg,    mesh0, [('outputnode.out_surf', 'surface2')]),
        (inputnode,  csv0, [('subject_id', 'subject_id')]),
        (mesh0,      csv0, [('distance', 'surf_dist')])
    ])

    mesh1 = pe.MapNode(ComputeMeshWarp(),
                       iterfield=['surface1', 'surface2'],
                       name='FMBSurfDistance')
    csv1 = pe.Node(AddCSVRow(in_file=settings['out_csv']),
                   name="FMBAddRow")
    csv1.inputs.method = 'FMB'

    wf.connect([
        (st1,       mesh1, [('out_dis_set.surf', 'surface1')]),
        (wrpsurf,   mesh1, [('out_points', 'surface2')]),
        (inputnode,  csv1, [('subject_id', 'subject_id')]),
        (mesh1,      csv1, [('distance', 'surf_dist')])
    ])

    mapen = ev.map_energy()
    wf.connect([
        (inputnode, mapen, [('subject_id', 'inputnode.subject_id')]),
        (regseg,    mapen, [('outputnode.out_enh', 'inputnode.reference'),
                            ('outputnode.reg_msk', 'inputnode.in_mask')]),
        (st1,       mapen, [('out_dis_set.surf', 'inputnode.surfaces0'),
                            ('out_ref_set.surf', 'inputnode.surfaces1')])
    ])

    cmethod1 = sdc_t2b()
    selmask = pe.Node(niu.Select(index=[2]), name='SelectMask')
    export2 = pe.Node(ExportSlices(all_axis=True), name='ExportT2B')
    mesh2 = pe.MapNode(ComputeMeshWarp(),
                       iterfield=['surface1', 'surface2'],
                       name='T2BSurfDistance')
    csv2 = pe.Node(AddCSVRow(in_file=settings['out_csv']),
                   name="T2BAddRow")
    csv2.inputs.method = 'T2B'

    wf.connect([
        (st1,        selmask, [('out_ref_set.segs', 'inlist')]),
        (selmask,   cmethod1, [('out', 'inputnode.t2w_mask')]),
        (st1,       cmethod1, [
            ('out_dis_set.dwi', 'inputnode.in_dwi'),
            ('out_dis_set.dwi_mask', 'inputnode.dwi_mask'),
            ('out_ref_set.t2w_brain', 'inputnode.in_t2w'),
            ('out_ref_set.surf', 'inputnode.in_surf')]),
        (ds,        cmethod1, [('bval', 'inputnode.in_bval'),
                               ('mr_param', 'inputnode.in_param')]),
        (cmethod1,   export2, [('outputnode.out_surf', 'surfaces0')]),
        (st1,        export2, [('out_dis_set.surf', 'surfaces1')]),
        (dti,        export2, [('outputnode.fa', 'reference')]),
        (st1,          mesh2, [('out_dis_set.surf', 'surface1')]),
        (cmethod1,     mesh2, [('outputnode.out_surf', 'surface2')]),
        (inputnode,     csv2, [('subject_id', 'subject_id')]),
        (mesh2,         csv2, [('distance', 'surf_dist')])
    ])

    return wf
