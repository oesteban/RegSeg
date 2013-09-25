#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os.path as op
import os

import nipype.pipeline.engine as pe
import dmri as dmri
import nipype.interfaces.utility as niu
import nipype.workflows.dmri.fsl as wf


def evaluation_workflow( name="Evaluation" ):
    evaluation = pe.Workflow( name='ISBI2014_Evaluation' )
    
    inputnode = pe.Node( niu.IdentityInterface( fields=['in_dwi','in_dwi_mask',
                                                        'in_t2', 'in_t2_mask',
                                                        'in_bvec','in_bval',
                                                        'in_rois','wm_mask' ] ),
                         name='inputnode' )
    
    te_incr = 2.46e-3 # secs
    dwell_time = 0.77e-3 # secs
    
    dist1_wf = dmri.distortion_workflow( name='Distortion_Y' )
    dist1_wf.inputs.inputnode.te_incr = te_incr
    dist1_wf.inputs.inputnode.echospacing = dwell_time
    dist1_wf.inputs.inputnode.encoding_direction = 'y'

    dist2_wf = dmri.distortion_workflow( name='Distortion_Yrev' )
    dist2_wf.inputs.inputnode.te_incr = te_incr
    dist2_wf.inputs.inputnode.echospacing = dwell_time
    dist2_wf.inputs.inputnode.encoding_direction = 'y-'
            
            
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
        # Connect the distortion workflows
         (inputnode,      dist1_wf, [ ('in_dwi','inputnode.in_file'),('in_dwi_mask','inputnode.in_mask') ])
        ,(inputnode,      dist2_wf, [ ('in_dwi','inputnode.in_file'),('in_dwi_mask','inputnode.in_mask') ])
        # Correct with fieldmap and perform tractography
        ,(dist1_wf,        fmap_wf, [ ('outputnode.out_file','inputnode.in_file'),('outputnode.out_phdiff_map','inputnode.fieldmap_pha') ])
        ,(inputnode,       fmap_wf, [ ('in_dwi_mask','inputnode.in_mask'),('in_dwi_mask','inputnode.fieldmap_mag' ) ])
        # Correct with topup
        ,(dist1_wf,       topup_wf, [ ('outputnode.out_file','inputnode.in_file_dir') ])
        ,(dist2_wf,       topup_wf, [ ('outputnode.out_file','inputnode.in_file_rev') ])
        # Correct with registration
        ,(dist1_wf,       t2reg_wf, [ ('outputnode.out_file','inputnode.in_file'),('outputnode.out_mask','inputnode.in_mask_dwi') ])
        ,(inputnode,      t2reg_wf, [ ('in_t2','inputnode.in_t2'),('in_t2_mask','inputnode.in_mask_t2')])
        #
        # Tractographies
        #
        # Perform tractography on the original phantom
        ,(inputnode,      trk_gold, [ ('in_dwi','inputnode.in_dwi'),('wm_mask','inputnode.in_mask' ), ('in_rois','inputnode.roi_file'),
                                      ('in_bvec','inputnode.in_bvec'),('in_bval','inputnode.in_bval' ) ])
        # Perform tractography on the fieldmap corrected
        ,(inputnode,      trk_dist, [ ('wm_mask','inputnode.in_mask' ), ('in_rois','inputnode.roi_file'),
                                      ('in_bvec','inputnode.in_bvec'),('in_bval','inputnode.in_bval' ) ])
        ,(dist1_wf,       trk_dist, [ ('outputnode.out_file','inputnode.in_dwi' ) ])
        # Perform tractography on the fieldmap corrected
        ,(inputnode,      trk_fmap, [ ('wm_mask','inputnode.in_mask' ), ('in_rois','inputnode.roi_file'),
                                      ('in_bvec','inputnode.in_bvec'),('in_bval','inputnode.in_bval' ) ])
        ,(fmap_wf,        trk_fmap, [ ('outputnode.epi_corrected','inputnode.in_dwi' ) ])
        # Perform tractography on the topup corrected
        ,(inputnode,       trk_rev, [ ('wm_mask','inputnode.in_mask' ), ('in_rois','inputnode.roi_file'),
                                      ('in_bvec','inputnode.in_bvec'),('in_bval','inputnode.in_bval' ) ])
        ,(topup_wf,        trk_rev, [ ('outputnode.epi_corrected','inputnode.in_dwi' ) ])
        # Perform tractography on the t2-registration corrected
        ,(inputnode,      trk_t2re, [ ('wm_mask','inputnode.in_mask' ), ('in_rois','inputnode.roi_file'),
                                      ('in_bvec','inputnode.in_bvec'),('in_bval','inputnode.in_bval' ) ])
        ,(t2reg_wf,       trk_t2re, [ ('outputnode.epi_corrected','inputnode.in_dwi' ) ])
    ])
    return evaluation


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
    evaluation.base_dir = op.join( working_dir, 'evaluation' )
    evaluation.inputs.inputnode.in_dwi = op.join( model_dir, 'DWIS_dti-scheme_no-noise.nii.gz' )
    evaluation.inputs.inputnode.in_dwi_mask = op.join( model_dir, 'signal_mask.nii.gz' )
    evaluation.inputs.inputnode.in_t2 = op.join( model_dir, "t2_weighted.nii.gz" )
    evaluation.inputs.inputnode.in_t2_mask = op.join( model_dir, "tpms", "tpm_mask.nii.gz" )
    evaluation.inputs.inputnode.in_bvec = op.join( model_dir, 'dti-scheme.bvec' ) 
    evaluation.inputs.inputnode.in_bval = op.join( model_dir, 'dti-scheme.bval' )
    evaluation.inputs.inputnode.in_rois = op.join( model_dir,'rois', 'seeding_regions.nii.gz' )
    evaluation.inputs.inputnode.wm_mask = op.join( model_dir, 'wm_mask.nii.gz' )
    evaluation.run()
    #evaluation.write_graph()


