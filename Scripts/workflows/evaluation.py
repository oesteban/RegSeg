#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os.path as op
import os

import nipype.pipeline.engine as pe
import dmri as dmri
import nipype.interfaces.utility as niu
import nipype.algorithms.misc as misc
import nipype.interfaces.nipy as nipy
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants
import nipype.workflows.dmri.fsl as wf


def evaluation_workflow( name="Evaluation" ):
    evaluation = pe.Workflow( name=name )
    
    inputnode = pe.Node( niu.IdentityInterface( fields=['in_dwi','in_dwi_mask',
                                                        'in_t2', 'in_t2_mask',
                                                        'in_bvec','in_bval',
                                                        'in_rois','wm_mask',
                                                        'out_csv','in_tpms' ] ),
                         name='inputnode' )
   
    outputnode = pe.Node( niu.IdentityInterface( fields=['out_tpms_resampled'] ),
                          name='outputnode' )
 
    te_incr = 2.46e-3 # secs
    dwell_time = 0.77e-3 # secs
    intensity = 0.75
    sigma = 8.0
    
    dist1_wf = dmri.distortion_workflow( name='Distortion_Y' )
    dist1_wf.inputs.inputnode.te_incr = te_incr
    dist1_wf.inputs.inputnode.echospacing = dwell_time
    dist1_wf.inputs.inputnode.encoding_direction = 'y'
    dist1_wf.inputs.inputnode.intensity = intensity
    dist1_wf.inputs.inputnode.sigma = sigma


    resample= pe.MapNode( fs.MRIConvert(out_type='niigz'), iterfield=['in_file'], name='Resample')
    normalize_tpms = pe.Node( niu.Function( input_names=['in_files'], output_names=['out_files'], function=dmri.normalize ), name='Normalize' )

    sel1 = pe.Node( niu.Select(index=3 ), name='SelectWM' )
    sel2 = pe.Node( niu.Select(index=3 ), name='SelectWM_dist' )
    sel3 = pe.Node( niu.Select(index=3 ), name='SelectWM_fmap' )
    sel4 = pe.Node( niu.Select(index=3 ), name='SelectWM_topup' )
    sel5 = pe.Node( niu.Select(index=3 ), name='SelectWM_t2' )

    sim1 = compare_dwis( name="Sim_dist" )
    sim1.inputs.inputnode.name='dist' 
    sim2 = compare_dwis( name="Sim_fmap" )
    sim2.inputs.inputnode.name='fmap' 
    sim3 = compare_dwis( name="Sim_reve" )
    sim3.inputs.inputnode.name='reve' 
    sim4 = compare_dwis( name="Sim_t2re" )
    sim4.inputs.inputnode.name='t2re' 

    ove1 = overlap_tpms_fmap()
    ove2 = overlap_tpms_topup()
    ove3 = overlap_tpms_t2reg()


    dist2_wf = dmri.distortion_workflow( name='Distortion_Yrev' )
    dist2_wf.inputs.inputnode.te_incr = te_incr
    dist2_wf.inputs.inputnode.echospacing = dwell_time
    dist2_wf.inputs.inputnode.encoding_direction = 'y-'
    dist2_wf.inputs.inputnode.intensity = intensity
    dist2_wf.inputs.inputnode.sigma = sigma
            
    fmap_wf = wf.fieldmap_correction()
    fmap_wf.inputs.inputnode.te_diff = te_incr*1e3 # in ms.
    fmap_wf.inputs.inputnode.epi_echospacing = dwell_time*1e3 # in ms
    fmap_wf.inputs.inputnode.vsm_sigma = 2.0
    fmap_wf.inputs.inputnode.encoding_direction = 'y'
    
    t2reg_wf = dmri.t2_registration_correct()
    
    topup_wf = wf.topup_correction()
    topup_wf.inputs.inputnode.encoding_direction = 'y'
    topup_wf.inputs.inputnode.ref_num = 0
    
    trk_gold = dmri.dtk_tractography_workflow(name='Tract_Gold')
    trk_fmap = trk_gold.clone('Tract_Fmap' )
    trk_t2re = trk_gold.clone('Tract_T2Reg' )
    trk_rev  = trk_gold.clone('Tract_Rev' )
    trk_dist = trk_gold.clone('Tract_Dist' )
    
    evaluation.connect([
        # Prepare tpms
         (inputnode,       resample, [ ('in_tpms','in_file'),(('in_dwi',dmri.get_vox_size),'vox_size') ])
        ,(resample,  normalize_tpms, [ ('out_file','in_files')])
        ,(normalize_tpms,outputnode, [ ('out_files','out_tpms_resampled')])
        ,(normalize_tpms,      sel1, [ ('out_files','inlist')])
        # Connect the distortion workflows
        ,(inputnode,      dist1_wf, [ ('in_dwi','inputnode.in_file'),('in_dwi_mask','inputnode.in_mask') ])
        ,(normalize_tpms, dist1_wf, [ ('out_files','inputnode.in_tpms') ])
        ,(inputnode,      dist2_wf, [ ('in_dwi','inputnode.in_file'),('in_dwi_mask','inputnode.in_mask') ])
        ,(normalize_tpms, dist2_wf, [ ('out_files','inputnode.in_tpms') ])
        ,(dist1_wf,           sel2, [ ('outputnode.out_tpms','inlist')])
        # Similarity after distortion
        ,(dist1_wf,           sim1, [ ('outputnode.out_file','inputnode.in_test') ])
        ,(inputnode,          sim1, [ ('in_dwi', 'inputnode.in_ref'),('in_dwi_mask','inputnode.in_mask'),('out_csv','inputnode.out_csv') ])
        # Correct with fieldmap
        ,(dist1_wf,        fmap_wf, [ ('outputnode.out_file','inputnode.in_file'),('outputnode.out_phdiff_map','inputnode.fieldmap_pha') ])
        ,(inputnode,       fmap_wf, [ ('in_dwi_mask','inputnode.in_mask'),('in_dwi_mask','inputnode.fieldmap_mag' ) ])
        # Similarity after fieldmap
        ,(fmap_wf,            sim2, [ ('outputnode.epi_corrected','inputnode.in_test') ])
        ,(inputnode,          sim2, [ ('in_dwi', 'inputnode.in_ref'),('in_dwi_mask','inputnode.in_mask'),('out_csv','inputnode.out_csv') ])
        # Correct with topup
        ,(dist1_wf,       topup_wf, [ ('outputnode.out_file','inputnode.in_file_dir') ])
        ,(dist2_wf,       topup_wf, [ ('outputnode.out_file','inputnode.in_file_rev') ])
        # Similarity after topup
        ,(topup_wf,           sim3, [ ('outputnode.epi_corrected','inputnode.in_test') ])
        ,(inputnode,          sim3, [ ('in_dwi', 'inputnode.in_ref'),('in_dwi_mask','inputnode.in_mask'),('out_csv','inputnode.out_csv') ])
        # Correct with registration
        ,(dist1_wf,       t2reg_wf, [ ('outputnode.out_file','inputnode.in_file'),('outputnode.out_mask','inputnode.in_mask_dwi') ])
        ,(inputnode,      t2reg_wf, [ ('in_t2','inputnode.in_t2'),('in_t2_mask','inputnode.in_mask_t2')])
        # Similarity after t2 registration
        ,(t2reg_wf,           sim4, [ ('outputnode.epi_corrected','inputnode.in_test') ])
        ,(inputnode,          sim4, [ ('in_dwi', 'inputnode.in_ref'),('in_dwi_mask','inputnode.in_mask'),('out_csv','inputnode.out_csv') ])
        #
        # Overlaps
        #
        ,(inputnode,          ove1, [ ('out_csv','inputnode.out_csv') ])
        ,(normalize_tpms,     ove1, [ ('out_files','inputnode.in_ref') ])
        ,(fmap_wf,            ove1, [ ('outputnode.out_vsm','inputnode.in_vsm') ])
        ,(dist1_wf,           ove1, [ ('outputnode.out_tpms','inputnode.in_tpms' ) ])
        ,(ove1,               sel3, [ ('outputnode.out_tpms','inlist')])
        ,(inputnode,          ove2, [ ('out_csv','inputnode.out_csv') ])
        ,(normalize_tpms,     ove2, [ ('out_files','inputnode.in_ref') ])
        ,(topup_wf,           ove2, [ ('outputnode.out_topup','inputnode.in_topup'),('outputnode.out_enc_file','inputnode.in_enc_file') ])
        ,(dist1_wf,           ove2, [ ('outputnode.out_tpms','inputnode.in_tpms' ) ])
        ,(dist2_wf,           ove2, [ ('outputnode.out_tpms','inputnode.in_tpms_rev') ])
        ,(ove2,               sel4, [ ('outputnode.out_tpms','inlist')])
        ,(inputnode,          ove3, [ ('out_csv','inputnode.out_csv') ])
        ,(normalize_tpms,     ove3, [ ('out_files','inputnode.in_ref') ])
        ,(t2reg_wf,           ove3, [ ('outputnode.fwd_tfm','inputnode.in_tfm') ])
        ,(dist1_wf,           ove3, [ ('outputnode.out_tpms','inputnode.in_tpms' ) ])
        ,(ove3,               sel5, [ ('outputnode.out_tpms','inlist')])
        #
        # Tractographies
        #
        # Perform tractography on the original phantom
        ,(inputnode,      trk_gold, [ ('in_dwi','inputnode.in_dwi'),('in_rois','inputnode.roi_file'),
                                      ('in_bvec','inputnode.in_bvec'),('in_bval','inputnode.in_bval' ) ])
        ,(sel1,           trk_gold, [ ('out','inputnode.in_mask') ])
        # Perform tractography on the distorted phantom
        ,(inputnode,      trk_dist, [ ('in_rois','inputnode.roi_file'),
                                      ('in_bvec','inputnode.in_bvec'),('in_bval','inputnode.in_bval' ) ])
        ,(dist1_wf,       trk_dist, [ ('outputnode.out_file','inputnode.in_dwi' ) ])
        ,(sel2,           trk_dist, [ ('out','inputnode.in_mask') ])
        # Perform tractography on the fieldmap corrected
        ,(inputnode,      trk_fmap, [ ('in_rois','inputnode.roi_file'),
                                      ('in_bvec','inputnode.in_bvec'),('in_bval','inputnode.in_bval' ) ])
        ,(fmap_wf,        trk_fmap, [ ('outputnode.epi_corrected','inputnode.in_dwi' ) ])
        ,(sel3,           trk_fmap, [ ('out','inputnode.in_mask') ])
        # Perform tractography on the topup corrected
        ,(inputnode,       trk_rev, [ ('in_rois','inputnode.roi_file'),
                                      ('in_bvec','inputnode.in_bvec'),('in_bval','inputnode.in_bval' ) ])
        ,(topup_wf,        trk_rev, [ ('outputnode.epi_corrected','inputnode.in_dwi' ) ])
        ,(sel4,            trk_rev, [ ('out','inputnode.in_mask') ])
        # Perform tractography on the t2-registration corrected
        ,(inputnode,      trk_t2re, [ ('in_rois','inputnode.roi_file'),
                                      ('in_bvec','inputnode.in_bvec'),('in_bval','inputnode.in_bval' ) ])
        ,(t2reg_wf,       trk_t2re, [ ('outputnode.epi_corrected','inputnode.in_dwi' ) ])
        ,(sel5,           trk_t2re, [ ('out','inputnode.in_mask') ])
    ])
    return evaluation


def compare_dwis( name='CompareDWIs' ):
    wf = pe.Workflow( name=name )

    inputnode = pe.Node( niu.IdentityInterface( fields=['in_test','in_ref','in_mask','out_csv','name'] ), name='inputnode' )
    outputnode = pe.Node( niu.IdentityInterface( fields=['similarity'] ), name='outputnode' )
    similarity = pe.Node( nipy.Similarity(metric='crl1'), name='Similarity' )
    savecsv = pe.Node( niu.Function( input_names=['fname','row','rowname','feature'], output_names=[], function=_append_csv ), name='AddCSVrow' )
    savecsv.inputs.feature = 'Similarity'

    wf.connect([
         (inputnode,   similarity, [('in_test','volume2'),('in_ref','volume1'),('in_mask','mask1'),('in_mask','mask2') ] )
        ,(inputnode,      savecsv, [('out_csv','fname'),('name','rowname')])
        ,(similarity,  outputnode, [('similarity','similarity') ])
        ,(similarity,     savecsv, [('similarity','row') ])
    ])
    return wf


def overlap_tpms_fmap( name='OverlapTpmsFmap' ):
    wf = pe.Workflow( name=name )

    inputnode= pe.Node( niu.IdentityInterface( fields=['in_tpms','in_ref','in_vsm','out_csv' ] ), name='inputnode' )
    outputnode = pe.Node( niu.IdentityInterface( fields=['out_tpms'] ), name='outputnode' )
    
    fixTpms = pe.MapNode( fsl.FUGUE(icorr=False), iterfield=['in_file'], name='FixTpms' )
    norm_tpms = pe.Node( niu.Function( input_names=['in_files'], output_names=['out_files'], function=dmri.normalize ), name='Normalize' )
    jaccard = pe.Node( misc.FuzzyOverlap(weighting='volume'), name='Jaccard' )
    merge = pe.Node( niu.Merge(2), name='merge' )
    savecsv = pe.Node( niu.Function( input_names=['fname','row','rowname','feature'], output_names=[], function=_append_csv ), name='AddCSVrow' )
    savecsv.inputs.feature='Overlap'
    savecsv.inputs.rowname='fmap'


    wf.connect([
         (inputnode,      fixTpms, [ ('in_tpms','in_file'),('in_vsm', 'shift_in_file' ) ])
        ,(inputnode,      jaccard, [ ('in_ref','in_ref' ) ])
        ,(fixTpms,      norm_tpms, [ ('unwarped_file','in_files') ])
        ,(norm_tpms,      jaccard, [ ('out_files','in_tst' ) ])
        ,(norm_tpms,   outputnode, [ ('out_files','out_tpms' ) ])
        ,(inputnode,      savecsv, [ ('out_csv','fname') ])
        ,(jaccard,          merge, [ ('jaccard','in1'),('class_fji','in2' ) ])
        ,(merge,          savecsv, [ ('out','row') ])
    ])
    return wf


def overlap_tpms_topup( name='OverlapTpmsTopup' ):
    wf = pe.Workflow( name=name )

    inputnode= pe.Node( niu.IdentityInterface( fields=['in_tpms','in_tpms_rev','in_ref','in_topup','in_enc_file','out_csv' ] ), name='inputnode' )
    outputnode = pe.Node( niu.IdentityInterface( fields=['out_tpms'] ), name='outputnode' )
    merg_tpms = pe.Node( fsl.Merge(dimension='t'), name='MergeTpms' )
    merg_rev = pe.Node( fsl.Merge(dimension='t'), name='MergeTpmsRev' )

    combin = pe.Node( niu.Merge(2), name='MergeEnc' )    
    applytopup = pe.Node( fsl.ApplyTOPUP(in_index=[1,2] ), name='applytopup' )

    split = pe.Node( fsl.Split( dimension='t' ), name='SplitEnc' )
    
    norm_tpms = pe.Node( niu.Function( input_names=['in_files'], output_names=['out_files'], function=dmri.normalize ), name='Normalize' )
    jaccard = pe.Node( misc.FuzzyOverlap(weighting='volume'), name='Jaccard' )
    merge = pe.Node( niu.Merge(2), name='merge' )
    savecsv = pe.Node( niu.Function( input_names=['fname','row','rowname','feature'], output_names=[], function=_append_csv ), name='AddCSVrow' )
    savecsv.inputs.feature='Overlap'
    savecsv.inputs.rowname='reve'


    wf.connect([
         (inputnode,    merg_tpms, [ ('in_tpms','in_files') ])
        ,(inputnode,     merg_rev, [ ('in_tpms_rev','in_files' ) ])
        ,(merg_tpms,       combin, [ ('merged_file','in1' ) ])
        ,(merg_rev,        combin, [ ('merged_file','in2' ) ])
        ,(inputnode,   applytopup, [ ('in_topup','in_topup'),('in_enc_file','encoding_file')])
        ,(combin,      applytopup, [ ('out','in_files')] )
        ,(applytopup,       split, [ ('out_corrected','in_file') ])
        ,(inputnode,      jaccard, [ ('in_ref','in_ref' ) ])
        ,(split,        norm_tpms, [ ('out_files','in_files') ])
        ,(norm_tpms,      jaccard, [ ('out_files','in_tst' ) ])
        ,(norm_tpms,   outputnode, [ ('out_files','out_tpms' ) ])
        ,(inputnode,      savecsv, [ ('out_csv','fname') ])
        ,(jaccard,          merge, [ ('jaccard','in1'),('class_fji','in2' ) ])
        ,(merge,          savecsv, [ ('out','row') ])
    ])

    return wf

def overlap_tpms_t2reg( name='OverlapTpmsT2Reg' ):
    wf = pe.Workflow( name=name )

    inputnode= pe.Node( niu.IdentityInterface( fields=['in_tpms','in_ref','in_tfm','out_csv' ] ), name='inputnode' )
    outputnode = pe.Node( niu.IdentityInterface( fields=['out_tpms'] ), name='outputnode' )
    
    select = pe.Node( niu.Select(index=[0]), name='select_ref' )

    applytfm = pe.MapNode( ants.WarpImageMultiTransform(), iterfield=['input_image'], name='ApplyTfm' ) 

    norm_tpms = pe.Node( niu.Function( input_names=['in_files'], output_names=['out_files'], function=dmri.normalize ), name='Normalize' )
    jaccard = pe.Node( misc.FuzzyOverlap(weighting='volume'), name='Jaccard' )
    merge = pe.Node( niu.Merge(2), name='merge' )
    savecsv = pe.Node( niu.Function( input_names=['fname','row','rowname','feature'], output_names=[], function=_append_csv ), name='AddCSVrow' )
    savecsv.inputs.feature='Overlap'
    savecsv.inputs.rowname='t2reg'


    wf.connect([
         (inputnode,       select, [ ('in_ref', 'inlist') ])
        ,(select,        applytfm, [ ('out','reference_image') ])
        ,(inputnode,     applytfm, [ ('in_tpms', 'input_image'),('in_tfm','transformation_series') ])
        ,(inputnode,      jaccard, [ ('in_ref','in_ref' ) ])
        ,(applytfm,     norm_tpms, [ ('output_image','in_files') ])
        ,(norm_tpms,      jaccard, [ ('out_files','in_tst' ) ])
        ,(norm_tpms,   outputnode, [ ('out_files','out_tpms' ) ])
        ,(inputnode,      savecsv, [ ('out_csv','fname') ])
        ,(jaccard,          merge, [ ('jaccard','in1'),('class_fji','in2' ) ])
        ,(merge,          savecsv, [ ('out','row') ])
    ])
    return wf



def _append_csv( fname, row, rowname, feature ):
    import csv
    row = [ rowname ] + [ feature ] + row
    with open( fname, 'a' ) as f:
        w = csv.writer( f )
        w.writerow( row )

if __name__ == '__main__':
    root_dir='/home/oesteban/workspace/ACWE-Reg/'
    data_dir= op.join( root_dir, 'Data' )
    mname = 'DigitalPhantomISBI'
    model_dir = op.join( data_dir, mname )
    tpms_dir = op.join( model_dir, 'tpms' )
    
    working_dir= op.join( root_dir, 'temp', 'ISBI2014' )
    
    snr_list = [ 'no-noise', 'SNR-10', 'SNR-20', 'SNR-30' ]
    
    if not op.exists( working_dir ):
        os.makedirs( working_dir )

    
    evaluation = evaluation_workflow()
    evaluation.base_dir = working_dir


    sim_file = op.join( working_dir, evaluation.name, 'results.csv' )
    if op.exists( sim_file ):
        os.remove( sim_file )

    evaluation.inputs.inputnode.in_dwi = op.join( model_dir, 'DWIS_dti-scheme_no-noise.nii.gz' )
    evaluation.inputs.inputnode.in_dwi_mask = op.join( model_dir, 'signal_mask.nii.gz' )
    evaluation.inputs.inputnode.in_t2 = op.join( model_dir, "t2_weighted.nii.gz" )
    evaluation.inputs.inputnode.in_t2_mask = op.join( model_dir, "tpms", "tpm_mask.nii.gz" )
    evaluation.inputs.inputnode.in_bvec = op.join( model_dir, 'dti-scheme.bvec' ) 
    evaluation.inputs.inputnode.in_bval = op.join( model_dir, 'dti-scheme.bval' )
    evaluation.inputs.inputnode.in_rois = op.join( model_dir,'rois', 'seeding_regions.nii.gz' )
    evaluation.inputs.inputnode.wm_mask = op.join( model_dir, 'wm_mask.nii.gz' )
    evaluation.inputs.inputnode.in_tpms = [ op.join( model_dir, 'tpms', '%s_volume_fraction.nii.gz' % n ) for n in [ 'bg','csf','gm','wm' ] ]
    evaluation.inputs.inputnode.out_csv = sim_file
    evaluation.run()
    #evaluation.write_graph()


