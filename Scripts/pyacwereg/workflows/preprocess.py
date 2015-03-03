#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2015-03-02 14:48:33
# @Last Modified by:   Oscar Esteban
# @Last Modified time: 2015-03-03 11:40:25


def preprocess(name='Preprocessing'):
    """
    Preprocessing of the regseg evaluation workflow
    """
    from nipype.pipeline import engine as pe
    from nipype.interfaces import utility as niu
    from nipype.interfaces import io as nio
    from nipype.interfaces import fsl

    from pyacwereg.workflows.fieldmap import bmap_registration
    from pysdcev.workflows.smri import all_surfaces

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

    ds_fields = ['t1w', 't2w', 't1w_brain', 't2w_brain', 't1w_mask', 'dwi',
                 'dwi_mask', 'bval', 'bvec', 'aseg', 'mr_param']
    out_fields = ds_fields + ['warped_dwi', 'warped_msk', 'warped_surf',
                              'bmap_wrapped', 'surf']

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['subject_id', 'data_dir']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=out_fields),
                         name='outputnode')

    ds_tpl_args = {k: [['subject_id', [v]]] for k, v in fnames.iteritems()}
    ds = pe.Node(nio.DataGrabber(
        infields=['subject_id'], outfields=ds_tpl_args.keys(),
        sort_filelist=False, template='*'), name='sMRISource')
    ds.inputs.field_template = {k: 'subjects/%s/%s'
                                for k in ds_tpl_args.keys()}
    ds.inputs.template_args = ds_tpl_args

    bmap_prep = bmap_registration()
    surfs = all_surfaces()

    dwisplit = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    wdwi = warp_dwi()
    dwimerge = pe.Node(fsl.Merge(dimension='t', output_type='NIFTI'),
                       name='MergeDWIs')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,  ds,        [('subject_id', 'subject_id'),
                                 ('data_dir', 'base_directory')]),
        (ds,         bmap_prep, [('t1w_brain', 'inputnode.t1w_brain'),
                                 ('dwi_mask', 'inputnode.dwi_mask'),
                                 ('fmap_mag', 'inputnode.mag'),
                                 ('fmap_pha', 'inputnode.pha')]),
        (ds,             surfs, [('t1w_brain', 'inputnode.norm'),
                                 ('dwi_mask', 'inputnode.in_mask'),
                                 ('aseg', 'inputnode.aseg')]),
        (ds,          dwisplit, [('dwi', 'in_file')]),
        (dwisplit,        wdwi, [('out_files', 'inputnode.in_dwis')]),
        (ds,              wdwi, [('dwi_mask', 'inputnode.dwi_mask'),
                                 ('mr_param', 'inputnode.mr_param')]),
        (bmap_prep,       wdwi, [('outputnode.unwrapped', 'inputnode.bmap')]),
        (surfs,           wdwi, [('outputnode.out_surf', 'inputnode.surf')]),
        (ds,        outputnode, [(f, f) for f in ds_fields]),
        (wdwi,        dwimerge, [('outputnode.dwis', 'in_files')]),
        (dwimerge,  outputnode, [('merged_file', 'warped_dwi')]),
        (wdwi,      outputnode, [('outputnode.dwi_mask', 'warped_msk'),
                                 ('outputnode.surf', 'warped_surf')]),
        (bmap_prep, outputnode, [('outputnode.wrapped', 'bmap_wrapped')]),
        (surfs,     outputnode, [('outputnode.out_surf', 'surf')])
    ])
    return wf


def warp_dwi(name='WarpDWIWorkflow'):
    """
    Distorts the splitted dwis given at the input following the
    theoretically - correct warping that corresponds to the bmap input
    """
    from nipype.pipeline import engine as pe
    from nipype.interfaces import utility as niu
    from nipype.interfaces import io as nio
    from nipype.algorithms.mesh import WarpPoints

    from pyacwereg.workflows.registration import apply_dfm
    from pyacwereg.workflows.fieldmap import vsm_fmb

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_dwis', 'dwi_mask', 'surf', 'mr_param', 'bmap']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['dwis', 'dwi_mask', 'surf']), name='outputnode')

    params = pe.Node(nio.JSONFileGrabber(), name='DWIparams')
    vsm = vsm_fmb(phase_unwrap=False)

    warp_data = apply_dfm(name='WarpData')
    warp_surf = pe.MapNode(WarpPoints(), iterfield=['points'],
                           name='WarpSurfs')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,     params, [('mr_param', 'in_file')]),
        (inputnode,        vsm, [('bmap', 'inputnode.in_bmap'),
                                 ('dwi_mask', 'inputnode.in_mask')]),
        (inputnode,  warp_surf, [('surf', 'points')]),
        (inputnode,  warp_data, [('in_dwis', 'inputnode.in_files'),
                                 ('dwi_mask', 'inputnode.in_mask')]),
        (params,           vsm, [('delta_te', 'inputnode.delta_te'),
                                 ('echospacing', 'inputnode.echospacing'),
                                 ('epi_acc', 'inputnode.acc_factor'),
                                 ('enc_dir', 'inputnode.enc_dir')]),
        (vsm,        warp_data, [('outputnode.dfm_inv', 'inputnode.dfm'),
                                 ('outputnode.jac_inv', 'inputnode.jac')]),
        (vsm,        warp_surf, [('outputnode.dfm', 'warp')]),

        (warp_data, outputnode, [('outputnode.out_files', 'dwis'),
                                 ('outputnode.out_mask', 'dwi_mask')]),
        (warp_surf, outputnode, [('out_points', 'surf')])
    ])
    return wf
