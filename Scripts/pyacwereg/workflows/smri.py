#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: Oscar Esteban - code@oscaresteban.es
# @Date:   2014-03-05 15:08:55
# @Last Modified by:   oesteban
# @Last Modified time: 2014-03-28 20:17:31

import os
import os.path as op

import nipype.interfaces.io as nio              # Data i/o
import nipype.interfaces.utility as niu         # utility
import nipype.pipeline.engine as pe             # pipeline engine
import nipype.interfaces.fsl as fsl             # fsl
import nipype.interfaces.freesurfer as fs       # freesurfer
import nipype.interfaces.ants as ants           # ANTS
import nipype.pipeline.engine as pe
import pyacwereg.utils.freesurfer as myfs


def prepare_smri( name='Prepare_sMRI'):
    """ A nipype workflow to generate structural MRI datasets, including:
    A T2w image, co-registered to the original T1w using bbreg
    T1w and T2w, Co-registered, brain-extracted, inhomogeneity corrected
    Surfaces set comprehending: right and left hemispheres both pial and white, and CSF.
    """
    wf = pe.Workflow( name=name )

    def _fsdir( path ):
        import os.path as op
        return op.join( path, 'FREESURFER' )

    inputnode = pe.Node( niu.IdentityInterface( fields=[ 'subject_id', 'data_dir' ] ), name='inputnode' )
    outputnode = pe.Node( niu.IdentityInterface( fields=[ 'out_surfs', 'out_smri',
                                                           'out_smri_brain', 'out_tpms' ] ),
                                                 name='outputnode' )

    ds = pe.Node( nio.DataGrabber(infields=['subject_id'], outfields=['t1w','t2w'], sort_filelist=False), name='DataSource' )
    ds.inputs.template = '*'
    ds.inputs.field_template = dict( t1w='subjects/%s/T1*.nii.gz', t2w='subjects/%s/T2*.nii.gz' )
    ds.inputs.template_args = dict( t1w=[['subject_id']], t2w=[['subject_id']] )

    fs_src = pe.Node( myfs.FSFiles(), name='FSSource')

    bbreg = pe.Node( fs.BBRegister( init='header', contrast_type='t2', registered_file=True ), name='T2w_to_T1w')
    T1toRAS = pe.Node( fs.MRIConvert( out_type='niigz', out_orientation='RAS' ), name='T1toRAS' )
    T2toRAS = pe.Node( fs.MRIConvert( out_type='niigz', out_orientation='RAS' ), name='T2toRAS' )
    merge_mri = pe.Node( niu.Merge(2), name='merge_mri')

    csf = csf_surface()
    surfs = surfs_to_native()
    merge_srf = pe.Node( niu.Merge(2), name='merge_surfs')
    tovtk = pe.MapNode( fs.MRIsConvert( out_datatype='vtk'), name='toVTK', iterfield=['in_file'])
    fixvtk = pe.MapNode( niu.Function( input_names=['in_file','in_ref'], output_names=['out_file'], function=myfs.fixvtk), name='fixVTK', iterfield=['in_file'])

    tfm_norm = pe.Node( fs.ApplyVolTransform(reg_header=True), name='norm_to_T1')
    T1brainToRAS = pe.Node( fs.MRIConvert( out_type='niigz', out_orientation='RAS' ), name='T1brainToRAS' )
    t2msk = pe.Node( fs.ApplyMask(), name='T2_BET' )
    n4_t2 = pe.Node( ants.N4BiasFieldCorrection( dimension=3 ), name='Bias_T2' )
    merge_brain = pe.Node( niu.Merge(2), name='merge_brain')

    fast = pe.Node( fsl.FAST( number_classes=3, img_type=1 ), name='SegmentT1' )

    wf.connect([
         ( inputnode,              ds, [ ('subject_id','subject_id'), ('data_dir','base_directory')])
        ,( inputnode,          fs_src, [ ('subject_id','subject_id'),
                                         (('data_dir',_fsdir),'fs_subjects_dir')])
        ,( inputnode,             csf, [ ('subject_id', 'inputnode.subject_id'),
                                         (('data_dir',_fsdir),'inputnode.fs_subjects_dir') ])
        ,( inputnode,           surfs, [ ('subject_id', 'inputnode.subject_id'),
                                         (('data_dir',_fsdir),'inputnode.fs_subjects_dir') ])
        ,( inputnode,           bbreg, [ ('subject_id', 'subject_id'),
                                         (('data_dir',_fsdir),'subjects_dir') ] )
        ,( ds,                  bbreg, [ ('t2w', 'source_file') ] )
        ,( ds,                T1toRAS, [ ('t1w', 'in_file')])
        ,( fs_src,           tfm_norm, [ ('norm','source_file') ])
        ,( T1toRAS,          tfm_norm, [ ('out_file','target_file') ])
        ,( T1toRAS,             surfs, [ ('out_file','inputnode.in_native')])
        ,( bbreg,             T2toRAS, [ ('registered_file','in_file')])
        ,( T1toRAS,           T2toRAS, [ ('out_file','reslice_like') ])
        ,( T1toRAS,         merge_mri, [ ('out_file','in1')])
        ,( T2toRAS,         merge_mri, [ ('out_file','in2')])
        ,( tfm_norm,     T1brainToRAS, [ ('transformed_file', 'in_file' )])
        ,( T1brainToRAS,  merge_brain, [ ('out_file', 'in1' )])
        ,( T1brainToRAS,        t2msk, [ ('out_file', 'mask_file') ])
        ,( T2toRAS,             t2msk, [ ('out_file', 'in_file')])
        ,( t2msk,               n4_t2, [ ('out_file', 'input_image')])
        ,( n4_t2,         merge_brain, [ ('output_image', 'in2')])
        ,( csf,             merge_srf, [ ('outputnode.out_surf', 'in1')])
        ,( surfs,           merge_srf, [ ('outputnode.out_surfs', 'in2')])
        ,( merge_srf,           tovtk, [ ('out', 'in_file')])
        ,( T1toRAS,              fast, [ ('out_file', 'in_files' ) ] )
        ,( T1toRAS,            fixvtk, [ ('out_file','in_ref')])
        ,( tovtk,              fixvtk, [ ('converted','in_file')])
        ,( merge_mri,      outputnode, [ ('out', 'out_smri')])
        ,( merge_brain,    outputnode, [ ('out', 'out_smri_brain')])
        ,( fixvtk,         outputnode, [ ('out_file','out_surfs')])
        ,( fast,           outputnode, [ ('partial_volume_files', 'out_tpms' ) ] )
    ])
    return wf


def extract_surface( name="GenSurface", labels=None ):
    """ A nipype workflow for surface extraction from labels in a segmentation.
    References:

    https://github.com/nipy/nipype/issues/307
    https://mail.nmr.mgh.harvard.edu/pipermail//freesurfer/2011-November/021391.html
    http://brainder.org/2012/05/08/importing-freesurfer-subcortical-structures-into-blender/
    https://mail.nmr.mgh.harvard.edu/pipermail/freesurfer/2013-June/030586.html
    """

    if labels is None:
        raise RuntimeError( "labels should contain an array of ids")
    pipeline = pe.Workflow( name=name )
    inputnode = pe.Node(niu.IdentityInterface( fields=['subject_id', 'fs_subjects_dir' ] ), name='inputnode' )
    outputnode = pe.Node(niu.IdentityInterface( fields=['out_surf', 'out_binary'  ] ), name='outputnode' )
    fs_src = pe.Node( myfs.FSFiles(), name='FreesurferSource')
    toRAS = pe.Node( fs.MRIConvert( out_type='mgz', out_orientation="RAS" ), name='toRAS' )
    tfm_norm = pe.Node( fs.ApplyVolTransform(reg_header=True), name='norm_to_native')
    label2vol = pe.Node( fs.Label2Vol(), name='aseg_to_native')
    binarize = pe.Node( niu.Function( input_names=['in_file', 'labels' ], output_names=["out_file"],
                                      function=myfs.merge_labels ), name='combine' )
    binarize.inputs.labels = labels

    convertmgz = pe.Node( fs.MRIConvert( out_type='mgz' ), name='convertmgz' )
    pretess = pe.Node( niu.Function( input_names=['in_file', 'in_norm', 'label' ], output_names=['out_file'],
                                    function=myfs.MRIPreTess ), name='pretess' )
    pretess.inputs.label = 1
    tess = pe.Node( fs.MRITessellate(label_value=1), name='tess' )
    smooth = pe.Node( fs.SmoothTessellation(disable_estimates=True ), name='mris_smooth' )
    rename = pe.Node( niu.Rename(keep_ext=False,format_string='surf.native'), name='rename')

    def _default_labels( in_labels ):
        from nipype.interfaces.base import isdefined
        if not isdefined( in_labels) or len( in_labels ) == 0:
            in_labels = [ 4, 5, 43, 44, 14, 15, 72, 24 ]
        return in_labels

    pipeline.connect( [
                        ( inputnode,    fs_src, [( 'subject_id', 'subject_id'), ('fs_subjects_dir','fs_subjects_dir') ])
                       ,( fs_src,        toRAS, [ ('rawavg', 'in_file') ])
                       ,( fs_src,    label2vol, [ ('aseg','seg_file'), ('aseg', 'reg_header')] )
                       ,( fs_src,     tfm_norm, [ ('norm','source_file') ])
                       ,( toRAS,      tfm_norm, [ ('out_file','target_file') ])
                       ,( toRAS,     label2vol, [ ('out_file', 'template_file') ])
                       ,( label2vol,  binarize, [ ('vol_label_file','in_file' ) ])
                       ,( binarize, convertmgz, [ ('out_file','in_file')])
                       ,( convertmgz,  pretess, [ ('out_file', 'in_file' )])
                       ,( tfm_norm,    pretess, [ ('transformed_file', 'in_norm' )])
                       ,( pretess,        tess, [ ('out_file', 'in_file' ) ])
                       ,( tess,         smooth, [ ('surface', 'in_file' ) ])
                       ,( smooth,       rename, [ ('surface', 'in_file' ) ])
                       ,( rename,   outputnode, [ ('out_file', 'out_surf' ) ])
                       ,( binarize, outputnode, [ ('out_file', 'out_binary')])
                       ])

    return pipeline

def csf_surface( name="CSF_Surface" ):
    """ An extract_surface workflow with labels set to liquid areas in aseg segmentation  """
    return extract_surface( name=name, labels=[ 4, 5, 43, 44, 14, 15, 72, 24 ] )

def surfs_to_native( name='Surfaces_to_native' ):
    """ A nipype workflow to project a surface from Freesurfer to the T1's native space,
    in vtk format.
    """
    wf = pe.Workflow(name=name)

    inputnode = pe.Node( niu.IdentityInterface( fields=['in_native','fs_subjects_dir','subject_id']), name='inputnode')
    get_tfm = pe.Node( niu.Function( input_names=['in_file', 'fs_subjects_dir', 'subject_id'],
                                     output_names=['out_tfm'],
                                     function=myfs.gen_fs_transform), name='surftfm')
    applytfm = pe.MapNode( niu.Function( input_names=['in_surf','in_reg','in_target','fs_subjects_dir','subject_id'],
                                         output_names=['out_surf'], function=myfs.transform_surface ),
                                         iterfield=['in_surf'], name='TransformSurface' )
    applytfm.inputs.in_surf = [ 'lh.white', 'rh.white','lh.pial', 'rh.pial' ]
    outputnode = pe.Node( niu.IdentityInterface( fields=['out_surfs'] ), name='outputnode' )

    wf.connect([
             ( inputnode,  get_tfm, [ ('in_native', 'in_file'), ('fs_subjects_dir', 'fs_subjects_dir'), ('subject_id', 'subject_id' )])
            ,( inputnode, applytfm, [ ('in_native', 'in_target'), ('fs_subjects_dir', 'fs_subjects_dir'), ('subject_id', 'subject_id' ) ])
            ,( get_tfm,   applytfm, [ ('out_tfm','in_reg') ])
            ,( applytfm, outputnode, [ ('out_surf','out_surfs') ])
        ])

    return wf
