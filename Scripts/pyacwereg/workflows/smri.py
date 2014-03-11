# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os
import os.path as op

import nipype.interfaces.io as nio              # Data i/o
import nipype.interfaces.utility as niu         # utility
import nipype.pipeline.engine as pe             # pypeline engine
import nipype.interfaces.fsl as fsl             # fsl
import nipype.interfaces.freesurfer as fs       # freesurfer
import nipype.interfaces.ants as ants           # ANTS
import nipype.pipeline.engine as pe
import pyacwereg.utils.freesurfer as myfs


def prepare_smri( name='Prepare_sMRI'):
    wf = pe.Workflow( name=name )

    def _fsdir( path ):
        import os.path as op
        return op.join( path, 'FREESURFER' )

    inputnode = pe.Node( niu.IdentityInterface( fields=[ 'subject_id', 'data_dir' ] ), name='inputnode' )
    outputnode = pe.Node( niu.IdentityInterface( fields=[ 'out_surfs', 'out_smri', 'out_smri_brain' ] ), name='outputnode' )

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
    n4 = pe.Node( ants.N4BiasFieldCorrection( dimension=3 ), name='T2_Bias' )
    merge_brain = pe.Node( niu.Merge(2), name='merge_brain')

    wf.connect([
                     ( inputnode,              ds, [ ('subject_id','subject_id'), ('data_dir','base_directory')])
                    ,( inputnode,          fs_src, [ ('subject_id','subject_id'), (('data_dir',_fsdir),'fs_subjects_dir')]) 
                    ,( inputnode,             csf, [ ('subject_id', 'inputnode.subject_id'), (('data_dir',_fsdir),'inputnode.fs_subjects_dir') ])
                    ,( inputnode,           surfs, [ ('subject_id', 'inputnode.subject_id'), (('data_dir',_fsdir),'inputnode.fs_subjects_dir') ])
                    ,( inputnode,           bbreg, [ ('subject_id', 'subject_id'), (('data_dir',_fsdir),'subjects_dir') ] )
                    ,( ds,                  bbreg, [ ('t2w', 'source_file') ] )
                    ,( ds,                T1toRAS, [ ('t1w', 'in_file')])
                    ,( fs_src,           tfm_norm, [ ('norm','source_file') ])
                    ,( T1toRAS,          tfm_norm, [ ('out_file','target_file') ])
                    ,( T1toRAS,             surfs, [ ('out_file','inputnode.in_native')])
                    ,( bbreg,             T2toRAS, [ ('registered_file','in_file')])
                    ,( T1toRAS,         merge_mri, [ ('out_file','in1')])
                    ,( T2toRAS,         merge_mri, [ ('out_file','in2')])
                    ,( tfm_norm,     T1brainToRAS, [ ('transformed_file', 'in_file' )])
                    ,( T1brainToRAS,  merge_brain, [ ('out_file', 'in1' )])
                    ,( T1brainToRAS,        t2msk, [ ('out_file', 'mask_file') ])
                    ,( T2toRAS,             t2msk, [ ('out_file', 'in_file')])
                    ,( t2msk,                  n4, [ ('out_file', 'input_image')])
                    ,( n4,            merge_brain, [ ('output_image', 'in2')])
                    ,( csf,             merge_srf, [ ('outputnode.out_surf', 'in1')])
                    ,( surfs,           merge_srf, [ ('outputnode.out_surfs', 'in2')])
                    ,( merge_srf,           tovtk, [ ('out', 'in_file')])
                    ,( T1toRAS,            fixvtk, [ ('out_file','in_ref')])
                    ,( tovtk,              fixvtk, [ ('converted','in_file')])
                    ,( merge_mri,      outputnode, [ ('out', 'out_smri')])
                    ,( merge_brain,    outputnode, [ ('out', 'out_smri_brain')])
                    ,( fixvtk,         outputnode, [ ('out_file','out_surfs')])
        ])
    return wf


def fieldmap_preparation( name="Fmap_prepare" ):
    workflow = pe.Workflow(name=name)
    
    # Setup i/o
    inputnode = pe.Node( niu.IdentityInterface( fields=[ 'subject_id', 'in_fmap_mag', 'in_fmap_pha', 'in_t1w_brain', 'in_surfs', 'in_t2w', 'fs_subjects_dir' ]), name='inputnode' )
    outputnode = pe.Node( niu.IdentityInterface(fields=[ 'out_fmap_mag', 'out_fmap_pha', 'out_smri', 'out_tpms', 'out_mask', 'out_surfs' ]), name='outputnode' )
    
    # Setup initial nodes
    fslroi = pe.Node( fsl.ExtractROI(t_min=0, t_size=1), name='GetFirst' )
    bet = pe.Node( fsl.BET( frac=0.4 ), name='BrainExtraction' )
    bbreg = pe.Node( fs.BBRegister( init='header', contrast_type='t2', registered_file=True ), name='T2w_to_T1w')
    applymsk = pe.Node( fs.ApplyMask(), name='MaskBrain_T2w' )
    n4 = pe.Node( ants.N4BiasFieldCorrection( dimension=3 ), name='Bias' )
    #windorise = pe.Node( fsl.ImageStats(op_string='-p 98'), name='Windorise' )
    tonifti0 = pe.Node( fs.MRIConvert(out_type="niigz", out_orientation="RAS" ), name='To_Nifti_0' )
    tonifti1 = pe.Node( fs.MRIConvert(out_type="niigz", out_orientation="RAS" ), name='To_Nifti_1' )
    tonifti2 = pe.Node( fs.MRIConvert(out_type="niigz", out_orientation="RAS" ), name='To_Nifti_2' )
    tonifti3 = pe.Node( fs.MRIConvert(out_type="niigz", out_orientation="RAS" ), name='To_Nifti_3' )
    
    # Setup ANTS and registration
    def _aslist( tname ):
        import numpy as np
        return np.atleast_1d( tname ).tolist()
    
    registration = pe.Node( ants.Registration(output_warped_image=True), name="FM_to_T1" )
    registration.inputs.transforms = ['Rigid','Affine','SyN'] #, 'SyN']
    registration.inputs.transform_parameters = [(2.0,),(1.0,),(0.75,4.0,2.0)] #,(0.2,1.0,1.0)]
    registration.inputs.number_of_iterations = [[50],[20],[15]] #,[10] ]
    registration.inputs.dimension = 3
    registration.inputs.metric = ['Mattes','CC', 'CC'] #,'CC']
    registration.inputs.metric_weight = [1.0]*3
    registration.inputs.radius_or_number_of_bins = [32,3,3] #,2]
    registration.inputs.sampling_strategy = ['Regular','Random','Random'] #,'Random']
    registration.inputs.sampling_percentage = [None,0.1,0.15] #,0.15]
    registration.inputs.convergence_threshold = [1.e-5,1.e-6,1.e-7] #,1.e-8]
    registration.inputs.convergence_window_size = [20,10,4] # ,3]
    registration.inputs.smoothing_sigmas = [[6.0],[4.0],[2.0]] #,[1.0]]
    registration.inputs.sigma_units = ['vox']*3
    registration.inputs.shrink_factors = [[6],[2],[2]] #,[1] ]
    registration.inputs.use_estimate_learning_rate_once = [True]*3
    registration.inputs.use_histogram_matching = [True]*3
    registration.inputs.initial_moving_transform_com = 0
    registration.inputs.collapse_output_transforms = True
    registration.inputs.winsorize_lower_quantile = 0.005
    registration.inputs.winsorize_upper_quantile = 0.975
    
    binarize = pe.Node( fs.Binarize( min=0.001 ), name='Binarize' )
    applyAnts = pe.Node( ants.ApplyTransforms(dimension=3,interpolation='BSpline' ), name='ApplyANTs' )
    msk_mag = pe.Node( fs.ApplyMask(), name='Brain_FMap_Mag' )
    bin_mag = pe.Node( fs.Binarize( min=0.001 ), name='OutMask' )
    # msk_pha = pe.Node( fs.ApplyMask(), name='Brain_FMap_Pha' )
    # smooth = pe.Node( fsl.SpatialFilter( operation='median', kernel_shape='sphere', kernel_size=4 ), name='SmoothPhaseMap' )
    fix_pha = pe.Node( niu.Function( input_names=['in_file'], output_names=['out_file'], function=check_range ), name='CheckPhaseMap' )
    n4_t1 = pe.Node( ants.N4BiasFieldCorrection( dimension=3 ), name='Bias_T1' )
    fast = pe.Node( fsl.FAST( number_classes=3, probability_maps=True, img_type=1 ), name='SegmentT1' )
    merge = pe.Node( niu.Merge(2), name='JoinNames')
    combine = pe.Node( fsl.Merge(dimension='t' ), name='Combine')
    
    workflow.connect([
                        # Connect inputs to nodes
                         ( inputnode,        tonifti0, [ ('in_fmap_mag', 'in_file')] )
                        ,( inputnode,        tonifti1, [ ('in_fmap_pha', 'in_file')] )
                        ,( inputnode,           bbreg, [ ('subject_id', 'subject_id'), ('in_t2w', 'source_file'),('fs_subjects_dir','subjects_dir') ] )
                        ,( inputnode,        tonifti3, [ ('in_t1w_brain', 'in_file' ) ] )
                        ,( inputnode,        applymsk, [ ('in_t1w_brain', 'mask_file') ] )
                        # Connections between nodes
                        ,( tonifti0,           fslroi, [ ('out_file', 'in_file')] )
                        ,( fslroi,                 n4, [ ('roi_file', 'input_image' ) ] )
                        ,( n4,                    bet, [ ('output_image','in_file' ) ] )
                        ,( bbreg,            applymsk, [ ('registered_file','in_file' ) ] )
                        ,( applymsk,         tonifti2, [ ('out_file','in_file') ])
                        # ANTs
                        ,( tonifti3,            n4_t1, [ ('out_file', 'input_image' ) ] )
                        ,( n4_t1,                fast, [ ('output_image', 'in_files' ) ] )
                        ,( n4_t1,        registration, [ ('output_image', 'fixed_image' ) ] )
                        ,( tonifti3,         binarize, [ ('out_file', 'in_file' ) ] )
                        ,( bet,          registration, [ ('out_file', 'moving_image' ) ] )
                        ,( binarize,     registration, [ ('binary_file', 'fixed_image_mask' ) ] )
                        ,( tonifti3,        applyAnts, [ ('out_file', 'reference_image' ) ] )
                        ,( tonifti1,        applyAnts, [ ('out_file', 'input_image' ) ] )
                        ,( registration,    applyAnts, [ ('forward_transforms','transforms'),('forward_invert_flags','invert_transform_flags') ])
                        ,( registration,      msk_mag, [ ('warped_image', 'in_file' ) ] )
                        ,( tonifti3,          msk_mag, [ ('out_file', 'mask_file' ) ] )
                        ,( msk_mag,           bin_mag, [ ('out_file', 'in_file') ])
                        #,( applyAnts,         msk_pha, [ ('output_image', 'in_file' ) ] )
                        #,( bin_mag,           msk_pha, [ ('binary_file', 'mask_file' ) ] )
                        #,( msk_pha,            smooth, [ ('out_file', 'in_file' ) ] )
                        ,( applyAnts,         fix_pha, [ ('output_image', 'in_file' ) ] )
                        ,( n4_t1,               merge, [ ('output_image', 'in1' ) ])
                        ,( tonifti2,            merge, [ ('out_file', 'in2') ])
                        ,( merge,             combine, [ ('out', 'in_files')])
                        # Connections to output
                        ,( combine,        outputnode, [ ('merged_file','out_smri') ])
                        ,( msk_mag,        outputnode, [ ('out_file', 'out_fmap_mag' ) ] )
                        ,( fix_pha,        outputnode, [ ('out_file', 'out_fmap_pha' ) ] )
                        ,( bin_mag,        outputnode, [ ('binary_file', 'out_mask' ) ] )
                        ,( fast,           outputnode, [ ('partial_volume_files', 'out_tpms' ) ] )
                     ])
    
    return workflow

def csf_surface( name="CSF_Surface" ):
    return extract_surface( name=name, labels=[ 4, 5, 43, 44, 14, 15, 72, 24 ] )

def extract_surface( name="GenSurface", labels=None ):
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

def surfs_to_native( name='Surfaces_to_native' ):
    wf = pe.Workflow(name=name)
    
    inputnode = pe.Node( niu.IdentityInterface( fields=['in_native','fs_subjects_dir','subject_id']), name='inputnode')
    get_tfm = pe.Node( niu.Function( input_names=['in_file', 'fs_subjects_dir', 'subject_id'],
                                     output_names=['out_tfm'],
                                     function=myfs.gen_fs_transform), name='surftfm')
    applytfm = pe.MapNode( niu.Function( input_names=['in_surf','in_reg','in_target','fs_subjects_dir','subject_id'],
                                         output_names=['out_surf'], function=myfs.transform_surface ), 
                                         iterfield=['in_surf'], name='TransformSurface' )
    applytfm.inputs.in_surf = [ 'lh.pial', 'rh.pial', 'lh.white', 'rh.white' ]
    outputnode = pe.Node( niu.IdentityInterface( fields=['out_surfs'] ), name='outputnode' )

    wf.connect([
             ( inputnode,  get_tfm, [ ('in_native', 'in_file'), ('fs_subjects_dir', 'fs_subjects_dir'), ('subject_id', 'subject_id' )])
            ,( inputnode, applytfm, [ ('in_native', 'in_target'), ('fs_subjects_dir', 'fs_subjects_dir'), ('subject_id', 'subject_id' ) ])
            ,( get_tfm,   applytfm, [ ('out_tfm','in_reg') ])
            ,( applytfm, outputnode, [ ('out_surf','out_surfs') ])
        ])

    return wf

# Helper functions

def check_range( in_file, out_file=None ):
    import nibabel as nb
    import os.path as op
    import numpy as np

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_checked.nii.gz' % fname)

    im = nb.load( in_file )
    imdata = im.get_data()

    imdata = imdata - imdata.min()
    norm = ( 8192 / imdata.max() )
    imdata = ( norm * imdata ) - 4096

    nii = nb.Nifti1Image( imdata, im.get_affine(), im.get_header() )
    nb.save( nii, out_file )
    return out_file