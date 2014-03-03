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



def smri_preparation( name="sMRI_prepare" ):
	workflow = pe.Workflow(name=name)

	# Setup i/o
	inputnode = pe.Node( niu.IdentityInterface( fields=[ 'subject_id', 'in_fmap_mag', 'in_fmap_pha', 'in_t1w_brain', 'in_t2w', 'fs_subjects_dir' ]), name='inputnode' )
	outputnode = pe.Node( niu.IdentityInterface(fields=[ 'out_fmap_mag', 'out_fmap_pha', 'out_t2w_brain' ]), name='outputnode' )

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
	registration.inputs.transforms = ['Rigid','Affine','SyN', 'SyN']
	registration.inputs.transform_parameters = [(2.0,),(1.0,),(1.0,4.0,6.0),(0.5,2.0,2.0)]
	registration.inputs.number_of_iterations = [[200],[100],[75],[50] ]
	registration.inputs.dimension = 3
	registration.inputs.metric = ['CC']*4
	registration.inputs.metric_weight = [1.0]*4
	registration.inputs.radius_or_number_of_bins = [2]*4
	registration.inputs.sampling_strategy = ['Regular','Random','Random','Random']
	registration.inputs.sampling_percentage = [None,0.2,0.15,0.15]
	registration.inputs.convergence_threshold = [1.e-5,1.e-7,1.e-8,1.e-8]
	registration.inputs.convergence_window_size = [20,10,5,5]
	registration.inputs.smoothing_sigmas = [[6.0],[4.0],[2.0],[1.0]]
	registration.inputs.sigma_units = ['vox']*4
	registration.inputs.shrink_factors = [[4],[2],[2],[1] ]
	registration.inputs.use_estimate_learning_rate_once = [True]*4
	registration.inputs.use_histogram_matching = [True]*4
	registration.inputs.initial_moving_transform_com = 0
	registration.inputs.collapse_output_transforms = True
	registration.inputs.winsorize_upper_quantile = 0.975

	binarize = pe.Node( fs.Binarize( min=0.001 ), name='Binarize' )
	applyAnts = pe.Node( ants.ApplyTransforms(dimension=3,interpolation='BSpline' ), name='ApplyANTs' )
	applymsk2 = pe.Node( fs.ApplyMask(), name='MaskBrain_FMap' )
	applymsk3 = pe.Node( fs.ApplyMask(), name='MaskBrain_FMap2' )

	workflow.connect([
	                  # Connect inputs to nodes
	                   ( inputnode, tonifti0, [ ('in_fmap_mag', 'in_file')] )
	                  ,( inputnode, tonifti1, [ ('in_fmap_pha', 'in_file')] )
	                  ,( inputnode, bbreg,    [ ('subject_id', 'subject_id'), ('in_t2w', 'source_file'),('fs_subjects_dir','subjects_dir') ] )
	                  ,( inputnode, tonifti3, [ ('in_t1w_brain', 'in_file' ) ] )
	                  ,( inputnode, applymsk, [ ('in_t1w_brain', 'mask_file') ] )
	                  # Connections between nodes
	                  ,( tonifti0, fslroi,    [ ('out_file', 'in_file')] )
	                  ,( fslroi,  n4,         [ ('roi_file', 'input_image' ) ] )
	                  ,( n4,            bet,  [ ('output_image','in_file' ) ] )
	                  ,( bbreg, applymsk,     [ ('registered_file','in_file' ) ] )
	                  # ANTs
	                  ,( tonifti3,     registration, [ ('out_file', 'fixed_image' ) ] )
	                  ,( tonifti3,         binarize, [ ('out_file', 'in_file' ) ] )
	                  ,( bet,          registration, [ ('out_file', 'moving_image' ) ] )
	                  ,( binarize,     registration, [ ('binary_file', 'fixed_image_mask' ) ] )
	                  ,( tonifti3,        applyAnts, [ ('out_file', 'reference_image' ) ] )
	                  ,( tonifti1,        applyAnts, [ ('out_file', 'input_image' ) ] )
	                  ,( registration,    applyAnts, [ ('forward_transforms','transforms'),('forward_invert_flags','invert_transform_flags') ])
	                  ,( applyAnts,       applymsk2, [ ('output_image', 'in_file' ) ] )
	                  ,( tonifti3,        applymsk2, [ ('out_file', 'mask_file' ) ] )
	                  ,( registration,    applymsk3, [ ('warped_image', 'in_file' ) ] )
	                  ,( tonifti3,        applymsk3, [ ('out_file', 'mask_file' ) ] )
	                  # Connections to output
	                  ,( applymsk,         tonifti2, [ ('out_file','in_file') ])
	                  ,( applymsk2,      outputnode, [ ('out_file', 'out_fmap_pha' ) ] )
	                  ,( tonifti2,       outputnode, [ ('out_file','out_t2w_brain') ])
	                  ,( applymsk3,      outputnode, [ ('out_file', 'out_fmap_mag' ) ] )
	                 ])

	return workflow

def distortion_workflow(name="synthetic_distortion"):
    pipeline = pe.Workflow(name=name)
    
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file', 'in_mask', 'te_incr', 'echospacing','encoding_direction','intensity','sigma', 'in_tpms' ]), name='inputnode' )
    phmap_siemens = pe.Node( niu.Function(input_names=['in_file','intensity','sigma'],output_names=['out_file'], function=generate_siemens_phmap ), name='gen_siemens_PhaseDiffMap' )
    prepare = pe.Node( fsl.PrepareFieldmap(nocheck=True), name='fsl_prepare_fieldmap' )
    vsm = pe.Node(fsl.FUGUE(save_shift=True), name="gen_VSM")
    dm = pe.Node( niu.Function(input_names=['in_file','in_mask'],output_names=['out_file'], function=demean ), name='demean' )
    applyWarp = fugue_all_workflow()
    vsm_fwd_mask = pe.Node(fsl.FUGUE(forward_warping=True), name='Fugue_WarpMask')
    binarize = pe.Node( fs.Binarize( min=0.00001 ), name="binarize" )
    warpimages = pe.MapNode( fsl.FUGUE(forward_warping=True,icorr=True), iterfield=['in_file'], name='WarpImages' )
    normalize_tpms = pe.Node( niu.Function( input_names=['in_files'], output_names=['out_files'], function=normalize ), name='Normalize' )
#    fixsignal = pe.Node( niu.Function( input_names=['in_reference','in_distorted','in_mask' ], output_names=['out_distorted'], function=sanitize ), name='FixSignal' )


    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file', 'out_vsm', 'out_mask', 'out_phdiff_map', 'out_tpms' ]), name='outputnode' )
    
    pipeline.connect([
                       (inputnode,phmap_siemens, [('in_mask', 'in_file'),('intensity','intensity'),('sigma','sigma') ] )
                      ,(phmap_siemens,  prepare, [('out_file','in_phase') ])
                      ,(inputnode,      prepare, [('in_mask','in_magnitude'),(('te_incr',_sec2ms),'delta_TE')])
                      ,(prepare,            vsm, [('out_fieldmap','phasemap_file')])
                      ,(inputnode,          vsm, [('in_mask','in_file'),('in_mask','mask_file'),('te_incr','asym_se_time'),('echospacing','dwell_time')])
                      ,(vsm,                 dm, [('shift_out_file', 'in_file')])
                      ,(inputnode,           dm, [('in_mask','in_mask')])
                      ,(inputnode, vsm_fwd_mask, [('in_mask','in_file'),('in_mask','mask_file')])
                      ,(dm,        vsm_fwd_mask, [('out_file','shift_in_file')])
                      ,(dm,           applyWarp, [('out_file', 'inputnode.in_vsm')])
                      ,(inputnode,    applyWarp, [('in_file', 'inputnode.in_file'),('in_mask', 'inputnode.in_mask'),('encoding_direction','inputnode.unwarp_direction') ])
                      ,(dm,          outputnode, [('out_file', 'out_vsm')])
                      ,(applyWarp,   outputnode, [('outputnode.out_file', 'out_file')])
                      #,(applyWarp,    fixsignal, [('outputnode.out_file', 'in_distorted')])
                      #,(inputnode,    fixsignal, [('in_file','in_reference')])
                      #,(vsm_fwd_mask, fixsignal, [('warped_file','in_mask')])
                      #,(fixsignal,   outputnode, [('out_distorted','out_file')])
                      ,(vsm_fwd_mask,  binarize, [('warped_file','in_file')])
                      ,(binarize,    outputnode, [('binary_file','out_mask')])
                      ,(phmap_siemens,outputnode,[('out_file','out_phdiff_map') ])
                      ,(inputnode,  warpimages, [('in_tpms','in_file'),('in_mask', 'mask_file') ])
                      ,(dm      ,   warpimages, [('out_file','shift_in_file') ])
                      ,(warpimages,  normalize_tpms, [('warped_file','in_files')])
                      ,(normalize_tpms,  outputnode, [('out_files','out_tpms')])
    ])
    
    return pipeline

#
# AUXILIARY FUGUE ALL WORKFLOW ---------------------------------------------------------------------------------
#

def fugue_all_workflow(name="Fugue_WarpDWIs"):
    pipeline = pe.Workflow(name=name)
    
    def _split_dwi(in_file):
        import nibabel as nib
        import os.path as op
        out_files = []
        frames = nib.funcs.four_to_three(nib.load(in_file))
        name, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            name, _ = op.splitext(name)
        for i, frame in enumerate(frames):
            out_file = op.abspath('./%s_%03d.nii.gz' % (name, i))
            nib.save(frame, out_file)
            out_files.append(out_file)
        return out_files
    
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file', 'in_mask', 'in_vsm','unwarp_direction' ]), name='inputnode' )
    dwi_split = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_files'], function=_split_dwi), name='split_DWI')
    vsm_fwd = pe.MapNode(fsl.FUGUE(forward_warping=True, icorr=False), iterfield=['in_file'], name='Fugue_Warp')
    dwi_merge = pe.Node(fsl.utils.Merge(dimension='t'), name='merge_DWI')
    
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file' ]), name='outputnode' )
    pipeline.connect([
                       (inputnode,  dwi_split, [('in_file', 'in_file')])
                      ,(dwi_split,    vsm_fwd, [('out_files', 'in_file')])
                      ,(inputnode,    vsm_fwd, [('in_mask', 'mask_file'),('in_vsm','shift_in_file'),('unwarp_direction','unwarp_direction') ])
                      ,(vsm_fwd,    dwi_merge, [('warped_file', 'in_files')])
                      ,(dwi_merge, outputnode, [('merged_file', 'out_file')])
                    ])    
    
    return pipeline




#
# HELPER FUNCTIONS ------------------------------------------------------------------------------------
#

def get_vox_size( in_file ):
    import nibabel as nib
    pixdim = nib.load( in_file ).get_header()['pixdim']
    return (pixdim[1],pixdim[2],pixdim[3])

def normalize( in_files ):
    import nibabel as nib
    import numpy as np
    import os.path as op

    out_files = []
    
    imgs = [ nib.load(fim) for fim in in_files ]
    img_data = np.array( [ im.get_data() for im in imgs ] ).astype( 'f32' )
    img_data[img_data>1.0] = 1.0
    img_data[img_data<0.0] = 0.0
    weights = np.sum( img_data, axis=0 )

    img_data[0][weights==0] = 1.0
    weights[weights==0] = 1.0
    

    nib.save( nib.Nifti1Image( weights, im.get_affine(), im.get_header() ), op.abspath('w.nii.gz' ) )
    
    for i,finname in enumerate( in_files ):
        fname,fext = op.splitext( op.basename( finname ) )
        if fext == '.gz':
            fname,fext2 = op.splitext( fname )
            fext = fext2 + fext
         
        out_file = op.abspath( fname+'_norm'+fext )
        out_files+= [ out_file ]
        data = img_data[i] / weights
        hdr = imgs[i].get_header()
        hdr['data_type']= 16
        hdr.set_data_dtype( 'float32' )   
        nib.save( nib.Nifti1Image( data, imgs[i].get_affine(), hdr ), out_file )
    return out_files


def sanitize( in_reference, in_distorted, in_mask, out_file=None ):
    import nibabel as nib
    import numpy as np
    import os.path as op
    import scipy.ndimage as ndimage

    inc = 7e-2
    msk_data = nib.load( in_mask ).get_data()
    msk_data[msk_data<(1.0-inc)] = 0
    msk_data[msk_data>(1.0+inc)] = 0
    msk_data[msk_data!=0] = 1.0
    msk_data = ndimage.binary_opening( msk_data ).astype( 'uint8' )
    im = nib.load( in_distorted )
    im_data = im.get_data()
    ref_data = nib.load( in_reference ).get_data()
        
    for i in range(0,np.shape(im_data)[-1]):
        im_data[:,:,:,i] = ref_data[:,:,:,i] * msk_data + im_data[:,:,:,i] * (1-msk_data) 
       
    if out_file is None:
        fname, fext = op.splitext(op.basename(in_distorted))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_fixed.nii.gz' % fname)
    nib.save( nib.Nifti1Image( im_data, im.get_affine(), im.get_header() ), out_file )
    return out_file


def generate_siemens_phmap( in_file, intensity, sigma, w1=-0.5, w2=1.0, out_file=None ):
    import nibabel as nib
    from scipy.ndimage import (binary_erosion,binary_dilation)
    from scipy.ndimage.filters import gaussian_filter
    import numpy as np
    import os.path as op
    import math

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_siemens.nii.gz' % fname)


    msk = nib.load( in_file )
    mskdata = msk.get_data()
    imshape = msk.get_shape()
    boundary = binary_erosion( mskdata ).astype(np.dtype('u1'))
    boundary = (( mskdata - boundary )).astype(float)
    
    x = ( int(0.65*imshape[0]),int(0.35*imshape[0]) )
    y = ( int(0.4*imshape[1]), int(0.6*imshape[1]) )
    z = ( int(0.4*imshape[2]), int( imshape[2] ) )

    slice1 = np.zeros( shape=imshape )
    slice1[x[0],y[0]:y[1],z[0]:z[1]]  = 1
    slice1 = binary_dilation( slice1 * boundary ).astype( np.float32 )

    slice2 = np.zeros( shape=imshape )
    slice2[x[1],y[0]:y[1],z[0]:z[1]]  = 1
    slice2 = binary_dilation( slice2 * boundary ).astype( np.float32 )
    distortionfront = intensity * ( slice1 * w1 + slice2 * w2 )

    phasefield = gaussian_filter( distortionfront, sigma=sigma )
   # phasefield = np.ma.array( phasefield, mask=1-mskdata )
   # median = np.ma.median( phasefield )
   # phasefield = phasefield - median

   # maxvalue = np.ma.max( phasefield.reshape(-1) )
   # if np.fabs( np.ma.min(phasefield.reshape(-1) ) )>maxvalue:
   #     maxvalue =np.fabs( np.ma.min(phasefield.reshape(-1) ) )

   # normalizer = 0.001 / maxvalue
    normalizer = 1.0

   # noisesrc = 2.0 * np.random.random_sample( size=np.shape(phasefield) ) - 1.0
   # noisesrc[mskdata>0] = 0.0
   # phasefield[mskdata==0] = 0.0
   # phasefield= phasefield * normalizer + noisesrc
    phasefield = phasefield * 2047.5

    hdr = msk.get_header()
    hdr['data_type'] = 16
    hdr['qform_code'] = 2
    hdr.set_data_dtype( np.float32 )
    nib.save( nib.Nifti1Image( phasefield.astype( np.float32 ), msk.get_affine(), hdr ), out_file ) 
    return out_file


def generate_phmap( in_file, out_file=None ):
    import nibabel as nib
    from scipy.ndimage import binary_erosion
    from scipy.ndimage.filters import gaussian_filter
    import numpy as np
    import os.path as op
    import math

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_ph_unwrap_rads.nii.gz' % fname)

    msk = nib.load( in_file )
    data = msk.get_data().astype(float)
    data = np.pi * data * (1.0/2048)
    hdr = msk.get_header()
    hdr.set_data_dtype( 16 )
    nib.save( nib.Nifti1Image( data, msk.get_affine(), hdr ), out_file ) 
    return out_file
    
def rad2radsec( in_file, in_mask, te_incr=2.46e-3, out_file=None ):
    import nibabel as nib
    import numpy as np
    import os.path as op
    
    im = nib.load( in_file )
    msk = nib.load( in_mask ).get_data()
    data = np.ma.array( im.get_data(), mask = 1-msk )
    data[msk==1] = data[msk==1] / ( te_incr )
    mval = np.ma.median( data )
    
    data1 = np.zeros( shape=im.get_shape() )
    data2 = np.zeros( shape=im.get_shape() )
    
    data1[msk==1] = data[msk==1] - mval
    
    data = np.rollaxis ( np.array( [ data1, data2 ] ), 0, 4)    
    hdr = im.get_header()
    hdr['dim'][0]=4
    hdr['dim'][4]=2
    
    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_ph_unwrap_radsec.nii.gz' % fname)
        
    nib.save( nib.Nifti1Image( data, im.get_affine(), hdr ), out_file )
    
    return out_file

def genmask(in_file, out_file=None):
    import numpy as np
    import nibabel as nib
    import os.path as op
    
    im = nib.load(in_file)
    
    data = np.ones( shape=im.get_shape() )
    
    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_mask.nii.gz' % fname)
        
    nib.save( nib.Nifti1Image( data, im.get_affine(), im.get_header() ), out_file )
    
    return out_file

def demean( in_file, in_mask, out_file=None ):
    import numpy as np
    import nibabel as nib
    import os.path as op
    
    im = nib.load( in_file )
    msk = nib.load( in_mask ).get_data()
    
    data = np.ma.array( im.get_data(), mask = 1-msk )
    mval = np.ma.median( data )
    data[msk==1] = data[msk==1] - mval
    
    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_demeaned.nii.gz' % fname)
        
    nib.save( nib.Nifti1Image( data, im.get_affine(), im.get_header() ), out_file )
    
    return out_file

def _ms2sec(val):
    return val*1e-3

def _sec2ms(val):
    return val*1e3

def _aslist( val ):
    return [ val ]

def _pickupWarp( val ):
    return val[0]


def _split_dwi(in_file):
    import nibabel as nib
    import os
    out_files = []
    frames = nib.funcs.four_to_three(nib.load(in_file))
    name, fext = os.path.splitext(os.path.basename(in_file))
    if fext == '.gz':
        name, _ = os.path.splitext(name)
    for i, frame in enumerate(frames):
        out_file = os.path.abspath('./%s_%03d.nii.gz' % (name, i))
        nib.save(frame, out_file)
        out_files.append(out_file)
    return out_files
