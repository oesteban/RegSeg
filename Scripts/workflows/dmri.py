# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os.path as op

import nipype.interfaces.io as nio              # Data i/o
import nipype.interfaces.utility as niu         # utility
import nipype.pipeline.engine as pe             # pypeline engine
import nipype.interfaces.fsl as fsl             # fsl
import nipype.interfaces.freesurfer as fs       # freesurfer
import nipype.interfaces.ants as ants           # ANTS


import nipype.interfaces.diffusion_toolkit as dtk
import nipype.interfaces.cmtk as cmtk
import nipype.interfaces.camino as camino
import nipype.interfaces.mrtrix as mrt
import nipype.interfaces.camino2trackvis as cam2trk
from tracks import ConnectivityMatrix

def t2_registration_correct( name="T2_Registration" ):
    pipeline = pe.Workflow( name=name )
    
    inputnode = pe.Node( niu.IdentityInterface( fields=['in_file',
                                                        'in_t2',
                                                        'in_mask_dwi',
                                                        'in_mask_t2']
                         ), name='inputnode' )

    outputnode = pe.Node( niu.IdentityInterface( fields=['epi_corrected','out_mask','fwd_tfm']),
                          name='outputnode' )

    get_b0 = pe.Node(fsl.ExtractROI( t_size=1, t_min=0 ), name="getB0")

    mask_t2 = pe.Node( fsl.ApplyMask(), name='maskT2' )
    mask_b0 = pe.Node( fsl.ApplyMask(), name='maskB0' )

    resample= pe.Node( fs.MRIConvert(out_type='niigz'),name='ResampleT2')

    reg = pe.Node( ants.Registration() , name='B0-to-T2' )

    reg.inputs.transforms = ['SyN']
    reg.inputs.transform_parameters = [(0.80, 0.5, 1.0) ]
    reg.inputs.number_of_iterations = [[200,100, 50]]
    reg.inputs.dimension = 3
    reg.inputs.write_composite_transform = True

    reg.inputs.metric = ['Mattes']
    reg.inputs.metric_weight = [1] # Default (value ignored currently by ANTs)
    reg.inputs.radius_or_number_of_bins = [64]
    reg.inputs.sampling_strategy = ['Random']
    reg.inputs.sampling_percentage = [0.1]
    reg.inputs.convergence_threshold = [1.e-9]
    reg.inputs.convergence_window_size = [20,10,6]
    reg.inputs.smoothing_sigmas = [[2,1,0]]
    reg.inputs.sigma_units = ['vox']
    reg.inputs.shrink_factors = [[3,2,1]]
    reg.inputs.use_estimate_learning_rate_once = [True]
    reg.inputs.use_histogram_matching = [True] # This is the default
#    reg.inputs.output_warped_image = 'output_warped_image.nii.gz'
#    reg.inputs.output_transform_prefix = "output_"
#    reg.inputs.collapse_output_transforms = False

    applytfm = pe.Node( ants.WarpTimeSeriesImageMultiTransform(), name='ApplyTfm' )
    #jacobian = pe.Node( ants.JacobianDeterminant(dimension=3,use_log=0), name='Jacobian' )
    #dwi_split = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_files'], function=_split_dwi), name='split_DWI')
    #applyjacobian = pe.MapNode( ants.MultiplyImages(dimension=3), iterfield=['first_input'], name='ApplyJacobian')
    #dwi_merge = pe.Node(fsl.Merge(dimension='t'), name='merge_DWI')


    applytfm_msk = pe.Node( ants.WarpImageMultiTransform(use_nearest=True), name='ApplyTfmMsk' )

    pipeline.connect( [
                         (inputnode,      get_b0, [ ('in_file','in_file') ])
                        ,(inputnode,    resample, [ ('in_t2', 'in_file'), ])
                        ,(get_b0,       resample, [(('roi_file',get_vox_size),'vox_size') ])
                        #,(inputnode,   dwi_split, [ ('in_file','in_files') ]
                        ,(inputnode,     mask_b0, [ ('in_mask_dwi', 'mask_file' ) ])
                        ,(inputnode,     mask_t2, [ ('in_mask_t2', 'mask_file' ) ])
                        ,(get_b0,        mask_b0, [ ('roi_file','in_file') ])
                        ,(inputnode,     mask_t2, [ ('in_t2','in_file') ])
                        ,(mask_t2,           reg, [ ( ('out_file',_aslist),'fixed_image') ])
                        ,(mask_b0,           reg, [ ( ('out_file',_aslist),'moving_image') ])
                        ,(inputnode,    applytfm, [ ('in_file','input_image') ])
                        ,(resample,     applytfm, [ ('out_file','reference_image') ])
                        ,(reg,          applytfm, [ ('forward_transforms','transformation_series') ])
                        ,(applytfm,   outputnode, [ ('output_image','epi_corrected' ) ])
                        #,(applytfm,    dwi_split, [ ('output_image','in_file') ])
                        #,(reg,          jacobian, [ (('forward_transforms',_pickupWarp),'warp_file') ])
                        #,(jacobian,applyjacobian, [ ('jacobian_image','second_input') ])
                        #,(dwi_split,applyjacobian,[ ('out_files','first_input')])
                        #,(applyjacobian,dwi_merge,[ ('output_product_image', 'in_files') ])
                        #,(dwi_merge,  outputnode, [ ('merged_file','epi_corrected') ])
                        ,(inputnode,applytfm_msk, [ ('in_mask_dwi','input_image') ])
                        ,(resample, applytfm_msk, [ ('out_file','reference_image') ])
                        ,(reg,      applytfm_msk, [ ('forward_transforms','transformation_series') ])
                        ,(applytfm_msk,outputnode,[ ('output_image','out_mask') ])
                        ,(reg,        outputnode, [ ('forward_transforms','fwd_tfm') ])
                      ])

    return pipeline



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

    pipeline.connect([
    ])

    
    return pipeline


def mrtrix_tractography_workflow( name='MRTrix_CSD_tracks', thres=0.5, seed='wm' ):
    pipeline = pe.Workflow(name=name)

    inputnode = pe.Node(  niu.IdentityInterface( 
                          fields=['in_dwi','in_bvec','in_bval','in_mask',
                                  'roi_file' ] ), 
                          name='inputnode' )

    outputnode = pe.Node( niu.IdentityInterface(
                          fields=['track_file','matrix_tracks',
                                  'connectome','avglen_map' ] ),
                          name='outputnode' )


    wmmask = pe.Node( fs.Binarize( min=thres ), name='BinWM' )
    fsl2mrtrix = pe.Node(mrt.FSL2MRTrix(),name='fsl2mrtrix')
    est_resp = pe.Node(mrt.EstimateResponseForSH(),name='est_resp')
    est_resp.inputs.maximum_harmonic_order = 6
    csdeconv = pe.Node(mrt.ConstrainedSphericalDeconvolution(),name='csdeconv')
    csdeconv.inputs.maximum_harmonic_order = 6
    track = pe.Node(mrt.SphericallyDeconvolutedStreamlineTrack(),name='CSDstreamtrack')
    track.inputs.desired_number_of_tracks = 150000
    tck2trk = pe.Node(mrt.MRTrix2TrackVis(),name='tck2trk')
    matrix = pe.Node( ConnectivityMatrix(), name='BuildMatrix' )

    pipeline.connect([
                         (inputnode,    fsl2mrtrix, [("in_bvec", "bvec_file"),('in_bval','bval_file') ])
                        ,(inputnode,        wmmask, [('in_mask','in_file')])
                        ,(inputnode,      est_resp, [("in_dwi","in_file")])
                        ,(fsl2mrtrix,     est_resp, [("encoding_file","encoding_file")])
                        ,(wmmask,         est_resp, [("binary_file","mask_image")])
                        ,(inputnode,      csdeconv, [("in_dwi","in_file")])
                        ,(wmmask,         csdeconv, [("binary_file","mask_image")])
                        ,(est_resp,       csdeconv, [("response","response_file")])
                        ,(fsl2mrtrix,     csdeconv, [("encoding_file","encoding_file")])
                        ,(csdeconv,          track, [("spherical_harmonics_image","in_file")])
                        ,(track,           tck2trk, [('tracked','in_file')])
                        ,(inputnode,       tck2trk, [('in_dwi','image_file')])
                        ,(tck2trk,      outputnode, [('out_file', 'track_file') ])
                        ,(tck2trk,          matrix, [('out_file','in_file')])
                        ,(inputnode,        matrix, [('roi_file','roi_file')])
                        ,(matrix,       outputnode, [('out_conmat','connectome'),('out_lenmat','avglen_map'),('out_tracks','matrix_tracks') ])
                     ])


    if seed=='wm':
        pipeline.connect([ (wmmask,          track, [("binary_file","seed_file")]) ])
    else:
        seedbin = pe.Node( fs.Binarize( min=thres ), name='BinSM' )
        pipeline.connect([ 
                             (inputnode,    seedbin, [('roi_file','in_file')])
                            ,(seedbin,        track, [("binary_file","seed_file")])
                        ])

    return pipeline




def dtk_tractography_workflow( name='DTK_tractography', seed='wm' ):
    pipeline = pe.Workflow(name=name)

    inputnode = pe.Node(  niu.IdentityInterface( 
                          fields=['in_dwi','in_bvec','in_bval','in_mask',
                                  'roi_file' ] ), 
                          name='inputnode' )

    outputnode = pe.Node( niu.IdentityInterface(
                          fields=['tensor','track_file','smoothed_track_file',
                                  'ADC','FA','connectome','avglen_map','filter_tracks' ] ),
                          name='outputnode' )

    dtifit = pe.Node(dtk.DTIRecon(output_type='nii.gz'),name='dtifit')
    dtk_tracker = pe.Node(dtk.DTITracker(mask1_threshold=[0.002,3.0],invert_x=False,invert_y=True,random_seed=10), name="dtk_tracker")

    smooth_trk = pe.Node(dtk.SplineFilter(step_length=0.5), name="smooth_trk")
    matrix = pe.Node( ConnectivityMatrix(), name='BuildMatrix' )

    # trk2camino = pe.Node(cam2trk.Trackvis2Camino(), name="trk2camino")
    # connectome = pe.Node(camino.Conmat(), name="Connectivity")
    
    pipeline.connect([
                         (inputnode,        dtifit, [('in_dwi','DWI'),('in_bvec','bvecs'),('in_bval','bvals')])
                        ,(inputnode,   dtk_tracker, [('in_mask','mask1_file') ])
                        ,(dtifit,      dtk_tracker, [('tensor','tensor_file') ])
                        ,(dtk_tracker,  smooth_trk, [('track_file', 'track_file')])
                        ,(dtifit,       outputnode, [('tensor', 'tensor'),('ADC','ADC'),('FA','FA') ])
                        ,(dtk_tracker,  outputnode, [('track_file', 'track_file') ])
                        ,(smooth_trk,   outputnode, [('smoothed_track_file','smoothed_track_file')])
                        ,(smooth_trk,       matrix, [('smoothed_track_file','in_file')])
                        ,(inputnode,        matrix, [('roi_file','roi_file')])
                        ,(matrix,       outputnode, [('out_conmat','connectome'),('out_lenmat','avglen_map'),('out_tracks','filter_tracks') ])
                        # ,(smooth_trk,   trk2camino, [('smoothed_track_file','in_file' ) ])
                        # ,(trk2camino,   connectome, [('camino','in_file')] )
                        # ,(inputnode,    connectome, [('roi_file','roi_file')])
                        # ,(connectome,   outputnode, [('out_file','connectome')])
                     ])

    if seed=='gm':
        mask_seed_th = [0.5, 54]
        pipeline.connect([
                         (inputnode,   dtk_tracker, [('roi_file','mask_seed') ])
                        ])
    else:
        mask_seed_th = [0.001,3.0]
        pipeline.connect([
                         (inputnode,   dtk_tracker, [('in_mask','mask_seed') ])
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
