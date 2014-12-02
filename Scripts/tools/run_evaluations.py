#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban - code@oscaresteban.es
# @Date:   2014-04-04 19:39:38
# @Last Modified by:   oesteban
# @Last Modified time: 2014-12-01 13:37:54

__author__ = "Oscar Esteban"
__copyright__ = "Copyright 2013, Biomedical Image Technologies (BIT), \
                 Universidad Politécnica de Madrid"
__credits__ = "Oscar Esteban"
__license__ = "FreeBSD"
__version__ = "0.1"
__maintainer__ = "Oscar Esteban"
__email__ = "code@oscaresteban.es"
__status__ = "Prototype"

try:
    from enthought.etsconfig.api import ETSConfig
    ETSConfig.toolkit = 'null'
except:
    pass
import os


def hcp_workflow(name='HCP_TMI2015', settings={}):
    """
    The pysdcev evaluation workflow for the human connectome project (HCP)
    """
    from nipype.pipeline import engine as pe
    from nipype.interfaces import utility as niu
    from nipype.interfaces import io as nio
    from nipype.interfaces import freesurfer as fs
    from nipype.algorithms.mesh import P2PDistance, WarpPoints
    from nipype.algorithms.misc import AddCSVRow
    from nipype.workflows.dmri.fsl.artifacts import sdc_fmb
    # from pysdcev.utils import all2RAS
    from pyacwereg import misc as acwregmisc
    from pyacwereg.interfaces.utility import ExportSlices
    from pyacwereg.workflows.registration import regseg_wf
    from pyacwereg.workflows import evaluation as ev
    from pysdcev.workflows.fieldmap import bmap_registration
    from pysdcev.workflows.warpings import process_vsm
    from pysdcev.workflows.smri import preprocess_t2, preprocess_dwi
    from pysdcev.workflows.tractography import mrtrix_dti
    from pysdcev.stages.stage1 import stage1

    print 'Subjects=' + str(settings['subject_id'])
    print 'Fieldmaps=' + str(settings['bmap_id'])

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['subject_id', 'data_dir', 'bmap_id']), name='inputnode')
    inputnode.inputs.data_dir = settings['data_dir']
    inputnode.iterables = [('subject_id', settings['subject_id']),
                           ('bmap_id', settings['bmap_id'])]

    fnames = dict(t1w='T1w_acpc_dc_restore.nii.gz',
                  t1w_brain='T1w_acpc_dc_restore_brain.nii.gz',
                  t2w='T2w_acpc_dc_restore.nii.gz',
                  t2w_brain='T2w_acpc_dc_restore_brain.nii.gz',
                  t1w_mask='brainmask_fs.nii.gz',
                  dwi='dwi.nii.gz',
                  dwi_mask='dwi_brainmask.nii.gz',
                  bval='bvals',
                  bvec='bvecs',
                  aseg='aparc+aseg.nii.gz')
    # surf='*.surf.gii')

    ds_tpl_args = {k: [['subject_id', [v]]] for k, v in fnames.iteritems()}

    ds = pe.Node(nio.DataGrabber(
        infields=['subject_id'], outfields=ds_tpl_args.keys(),
        sort_filelist=False, template='*'), name='sMRISource')
    ds.inputs.field_template = {k: 'subjects/%s/%s'
                                for k in ds_tpl_args.keys()}
    ds.inputs.template_args = ds_tpl_args

    # surfsort = pe.Node(niu.Function(
    #     function=acwregmisc.sort_surfs, input_names=['surfs'],
    #     output_names=['out']), name='SurfSorted')

    bmapnames = dict(mag='FM_mag.nii.gz',
                     pha='FM_pha.nii.gz',
                     param='parameters.txt')
    bm_tpl_args = {k: [['bmap_id', [v]]] for k, v in bmapnames.iteritems()}

    ds_bmap = pe.Node(nio.DataGrabber(
        infields=['bmap_id'], outfields=bm_tpl_args.keys(),
        sort_filelist=False, template='*'), name='FieldmapSource')
    ds_bmap.inputs.field_template = {k: 'fieldmaps/%s/%s'
                                     for k in bm_tpl_args.keys()}
    ds_bmap.inputs.template_args = bm_tpl_args

    rfield = [k for k, f in fnames.iteritems() if '.nii' in f]
    rparam = {k: dict(resample_type='nearest') for k in rfield
              if (('brain' in k) or ('aseg' in k))}

    # reorient = all2RAS(input_fields=rfield, input_param=rparam)
    bmap_prep = bmap_registration(factor=6.0)
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
        (inputnode,  ds_bmap,   [('bmap_id', 'bmap_id'),
                                 ('data_dir', 'base_directory')]),
        (ds_bmap,    bmap_prep, [('mag', 'inputnode.mag'),
                                 ('pha', 'inputnode.pha')]),
        (ds,         poly_msk,  [('aseg', 'in_file')]),
        (poly_msk,   pmregrid,  [('binary_file', 'in_file')]),
        (ds,         pmregrid,  [('dwi_mask', 'reslice_like')]),
        (ds,         bmap_prep, [('t1w_brain', 'inputnode.t1w_brain'),
                                 ('dwi_mask', 'inputnode.dwi_mask')]),
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
        # (ds,   surfsort, [('surf', 'surfs')]),
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
                          ('aseg', 'inputnode.parcellation')]),
        # (surfsort,  st1, [('out', 'inputnode.surf')]),
        (ds_bmap,   st1, [('param', 'inputnode.mr_params')]),
        (bmap_prep, st1, [
            ('outputnode.wrapped', 'inputnode.bmap_wrapped'),
            ('outputnode.unwrapped', 'inputnode.bmap_unwrapped')])
    ])

    dti = mrtrix_dti()
    rlmsk = pe.Node(fs.MRIConvert(), name='MaskReslice')
    mdti = pe.Node(niu.Merge(2), name='MergeDTI')
    msk = pe.Node(niu.Select(index=[2]), name='SelectHIresMask')

    regseg = regseg_wf()
    regseg.inputs.inputnode.iterations = [150, 100, 100]
    regseg.inputs.inputnode.descript_update = [30, 30, 30]
    regseg.inputs.inputnode.step_size = [0.1, 0.2, 0.3]
    regseg.inputs.inputnode.alpha = [0.0, 0.0, 0.0]
    regseg.inputs.inputnode.beta = [0.0, 0.0, 0.0]
    regseg.inputs.inputnode.convergence_energy = [True, True, True]
    regseg.inputs.inputnode.convergence_window = [8, 10, 15]
    regseg.inputs.inputnode.convergence_value = [1.0e-5, 1.0e-8, 1.0e-9]
    regseg.inputs.inputnode.f_smooth = [2.0, 0.5, None]
    regseg.inputs.inputnode.images_verbosity = 3
    regseg.inputs.inputnode.scales = [(0.0, 1.0, 0.0)] * 3
    regseg.inputs.inputnode.grid_spacing = [
        (20., 45., 10.), (10., 20., 10.), (10., 10., 10.)]

    wf.connect([
        (st1,   dti,    [('out_dis_set.dwi', 'inputnode.in_dwi'),
                         ('out_dis_set.dwi_mask', 'inputnode.in_mask')]),
        (ds,    dti,    [('bvec', 'inputnode.in_bvec'),
                         ('bval', 'inputnode.in_bval')]),
        (dti,   mdti,   [('outputnode.fa', 'in1'),
                         ('outputnode.md', 'in2')]),
        (mdti,  regseg, [('out', 'inputnode.in_fixed')]),
        (st1,   regseg, [('out_dis_set.tpms', 'inputnode.in_tpms'),
                         ('out_ref_set.surf', 'inputnode.in_surf')]),
        (st1,   msk,    [('out_dis_set.segs', 'inlist')]),
        (msk,   rlmsk,  [('out', 'in_file')]),
        (dti,   rlmsk,  [('outputnode.fa', 'reslice_like')]),
        (rlmsk, regseg, [('out_file', 'inputnode.in_mask')])
    ])

    cmethod0 = sdc_fmb(bmap_params=dict(delta_te=2.46e-3),
                       epi_params=dict(echospacing=0.77e-3,
                                       acc_factor=3,
                                       enc_dir='y-'))
    selbmap = pe.Node(niu.Split(splits=[1, 1], squeeze=True),
                      name='SelectBmap')
    dfm = process_vsm()
    dfm.inputs.inputnode.scaling = 1.0
    dfm.inputs.inputnode.enc_dir = 'y-'
    sunwarp = pe.MapNode(WarpPoints(), iterfield=['points'],
                         name='UnwarpSurfs')

    wf.connect([
        (st1,       cmethod0, [('out_dis_set.dwi', 'inputnode.in_file'),
                               ('out_dis_set.dwi_mask', 'inputnode.in_mask')]),
        (ds,        cmethod0, [('bval', 'inputnode.in_bval')]),
        (bmap_prep,  selbmap, [('outputnode.wrapped', 'inlist')]),
        (selbmap,   cmethod0, [('out1', 'inputnode.bmap_mag'),
                               ('out2', 'inputnode.bmap_pha')]),
        (cmethod0,       dfm, [('outputnode.out_vsm', 'inputnode.vsm')]),
        (st1,            dfm, [
            ('out_dis_set.dwi_mask', 'inputnode.reference')]),
        (dfm,        sunwarp, [('outputnode.inv_dfm', 'warp')]),
        (st1,        sunwarp, [('out_dis_set.surf', 'points')])
    ])

    export0 = pe.Node(ExportSlices(all_axis=True), name='ExportREGSEG')
    export1 = pe.Node(ExportSlices(all_axis=True), name='ExportFMB')

    wf.connect([
        (regseg,   export0, [('outputnode.out_surf', 'surfaces0')]),
        (st1,      export0, [('out_dis_set.surf', 'surfaces1')]),
        (dti,      export0, [('outputnode.fa', 'reference')]),
        (sunwarp,  export1, [('out_points', 'surfaces0')]),
        (st1,      export1, [('out_ref_set.surf', 'surfaces1')]),
        (rdti,     export1, [('outputnode.fa', 'reference')])
    ])

#    mesh0 = pe.MapNode(P2PDistance(weighting='surface'),
#                       iterfield=['surface1', 'surface2'],
#                       name='REGSEGSurfDistance')
#    csv0 = pe.Node(AddCSVRow(in_file=settings['out_csv']),
#                   name="REGSEGAddRow")
#    csv0.inputs.method = 'REGSEG'
#
#    wf.connect([
#        (st1,       mesh0, [('out_dis_set.surf', 'surface1')]),
#        (regseg,    mesh0, [('outputnode.out_surf', 'surface2')]),
#        (inputnode,  csv0, [('subject_id', 'subject_id')]),
#        (mesh0,      csv0, [('distance', 'surf_dist')])
#    ])
#
#    mesh1 = pe.MapNode(P2PDistance(weighting='surface'),
#                       iterfield=['surface1', 'surface2'],
#                       name='FMBSurfDistance')
#    csv1 = pe.Node(AddCSVRow(in_file=settings['out_csv']),
#                   name="FMBAddRow")
#    csv1.inputs.method = 'FMB'
#
#    wf.connect([
#        (st1,       mesh1, [('out_ref_set.surf', 'surface1')]),
#        (sunwarp,   mesh1, [('out_points', 'surface2')]),
#        (inputnode,  csv1, [('subject_id', 'subject_id')]),
#        (mesh1,      csv1, [('distance', 'surf_dist')])
#    ])

    return wf

if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    import os.path as op
    from glob import glob

    parser = ArgumentParser(description='TMI2015 - Experiment on HCP data',
                            formatter_class=RawTextHelpFormatter)

    g_input = parser.add_argument_group('Input')
    g_input.add_argument('-S', '--subjects_dir', action='store',
                         default=os.getenv('NEURO_DATA_HOME', '..'),
                         help='directory where subjects should be found')
    g_input.add_argument('-s', '--subject', action='store', default='*',
                         help='subject id or pattern')
    g_input.add_argument('-w', '--work_dir', action='store',
                         default=os.getcwd(),
                         help='directory to store intermediate results')
    g_input.add_argument('-f', '--fieldmap_id', action='store',
                         default='*', help='fieldmap id or pattern')
    # g_input.add_argument('-t', '--tstep', action='store',
    #                      default='*', help='subject id or pattern')
    g_input.add_argument('-N', '--name', action='store',
                         default='TMI2015_EXP2',
                         help=('default workflow name, '
                               'it will create a new folder'))

    g_output = parser.add_argument_group('Output')
    g_output.add_argument('--out_csv', action='store',
                          help=('default output csv file'))

    g_options = parser.add_argument_group('Settings')
    g_options.add_argument('-T', '--type', action='store',
                           choices=['fmb', 'peb', 'fsl'], default='fmb',
                           help='select SDC workflow type')

    opts = parser.parse_args()

    if not op.exists(opts.work_dir):
        os.makedirs(opts.work_dir)

    settings = {}
    settings['work_dir'] = opts.work_dir
    settings['data_dir'] = op.abspath(opts.subjects_dir)
    settings['bmap_id'] = [
        op.basename(f) for f in glob(op.join(opts.subjects_dir, 'fieldmaps',
                                             opts.fieldmap_id))]

    subjects = [
        op.basename(sub) for sub in glob(op.join(opts.subjects_dir, 'subjects',
                                                 opts.subject))]
    settings['subject_id'] = subjects

    if opts.out_csv is None:
        settings['out_csv'] = op.join(opts.work_dir, opts.name, 'results.csv')
    else:
        settings['out_csv'] = opts.out_csv

    wf = hcp_workflow(name=opts.name, settings=settings)
    wf.base_dir = settings['work_dir']
    wf.write_graph(format='pdf')
    wf.run()
