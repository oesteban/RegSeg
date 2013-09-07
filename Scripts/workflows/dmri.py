# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os.path as op

import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as niu      # utility
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.fsl as fsl          # fsl
import nipype.interfaces.freesurfer as fs    # freesurfer
import nipype.interfaces.diffusion_toolkit as dtk
import nipype.interfaces.mrtrix as mrtrix    #<---- The important new part!
import nipype.interfaces.camino as camino

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
    
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file', 'in_mask', 'in_vsm' ]), name='inputnode' )
    dwi_split = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_files'], function=_split_dwi), name='split_DWI')
    vsm_fwd = pe.MapNode(fsl.FUGUE(forward_warping=True), iterfield=['in_file'], name='Fugue_Warp')
    dwi_merge = pe.Node(fsl.utils.Merge(dimension='t'), name='merge_DWI')
    
    outputnode = pe.Node(interface=niu.IdentityInterface(fields=['out_file' ]), name='outputnode' )
    pipeline.connect([
                       (inputnode,  dwi_split, [('in_file', 'in_file')])
                      ,(dwi_split,    vsm_fwd, [('out_files', 'in_file')])
                      ,(inputnode,    vsm_fwd, [('in_mask', 'mask_file'),('in_vsm','shift_in_file')])
                      ,(vsm_fwd,    dwi_merge, [('warped_file', 'in_files')])
                      ,(dwi_merge, outputnode, [('merged_file', 'out_file')])
                    ])    
    
    return pipeline

def distortion_workflow(name="synthetic_distortion"):
    pipeline = pe.Workflow(name=name)
    
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file', 'in_mask', 'te_incr', 'echospacing' ]), name='inputnode' )
    phmap = pe.Node( niu.Function(input_names=['in_file'],output_names=['out_file'], function=generate_phmap ), name='gen_PhaseDiffMap' )
    maskgen = pe.Node( niu.Function(input_names=['in_file'],output_names=['out_file'], function=genmask ), name='gen_Mask' )
    prelude = pe.Node(fsl.PRELUDE(process3d=True), name='PhaseUnwrap')
    prepare = pe.Node(niu.Function(input_names=['in_file','in_mask','te_incr'], output_names=['out_file'], function=rad2radsec ), name='PhaseRadSec' )
    vsm = pe.Node(fsl.FUGUE(save_shift=True), name="gen_VSM")
    dm = pe.Node( niu.Function(input_names=['in_file','in_mask'],output_names=['out_file'], function=demean ), name='demean' )
    applyWarp = fugue_all_workflow()
    vsm_fwd_mask = pe.Node(fsl.FUGUE(forward_warping=True), name='Fugue_WarpMask')
    binarize = pe.Node( fs.Binarize( min=0.00001 ), name="binarize" )
    
    outputnode = pe.Node(interface=niu.IdentityInterface(fields=['out_file', 'out_vsm', 'out_mask' ]), name='outputnode' )
    
    
    pipeline.connect([
                       (inputnode,      phmap, [('in_mask', 'in_file')] )
                      ,(inputnode,    maskgen, [('in_mask', 'in_file')] )
                      ,(inputnode,    prelude, [('in_mask','magnitude_file')] )
                      ,(phmap,        prelude, [('out_file', 'phase_file')] )
                      ,(prelude,      prepare, [('unwrapped_phase_file','in_file')])
                      ,(inputnode,    prepare, [('in_mask','in_mask'),('te_incr','te_incr')])
                      ,(maskgen,          vsm, [('out_file','mask_file')])
                      ,(prepare,          vsm, [('out_file','phasemap_file')])
                      ,(inputnode,        vsm, [('in_mask','in_file'),('te_incr','asym_se_time'),('echospacing','dwell_time')])
                      ,(vsm,               dm, [('shift_out_file', 'in_file')])
                      ,(inputnode,         dm, [('in_mask','in_mask')])
                      ,(inputnode,vsm_fwd_mask,[('in_mask','in_file'),('in_mask','mask_file')])
                      ,(dm,      vsm_fwd_mask, [('out_file','shift_in_file')])
                      ,(dm,         applyWarp, [('out_file', 'inputnode.in_vsm')])
                      ,(inputnode,  applyWarp, [('in_file', 'inputnode.in_file'),('in_mask', 'inputnode.in_mask')])
                      ,(dm,        outputnode, [('out_file', 'out_vsm')])
                      ,(applyWarp, outputnode, [('outputnode.out_file', 'out_file')])
                      ,(vsm_fwd_mask,binarize, [('warped_file','in_file')])
                      ,(binarize,  outputnode, [('binary_file','out_mask')])
                      ])
    
    return pipeline

def fsl_fitting_workflow(name="DTIFit", out_dir=None ):
    pipeline = pe.Workflow(name=name)

    if out_dir is None:
        out_dir = op.abspath( './' )
    
    inputnode = pe.Node( niu.IdentityInterface( fields=['in_file', 'in_mask', 
                'in_bvec','in_bval']), name='inputnode' )
    
    fsmask = pe.Node( fs.ApplyMask(), name='MaskRawData' )
    dtifit = pe.Node( fsl.DTIFit(), name='TensorFitting' )

    datasink = pe.Node( nio.DataSink( base_directory=out_dir,
                        parameterization=False, container='results'),
                        name='sinker')
    # connect nodes
    pipeline.connect([
                       (inputnode, fsmask, [('in_file','in_file'),('in_mask','mask_file') ])
                      ,(inputnode, dtifit, [('in_file','dwi'),('in_mask','mask'),('in_bvec','bvec'),('in_bval','bval') ])
                      ,(dtifit,datasink, [( 'FA','@FA'), ('MD','@ADC') ] )
                     ])
    return pipeline


#
# HELPER FUNCTIONS ------------------------------------------------------------------------------------
#

def generate_phmap( in_file, out_file=None ):
    import nibabel as nib
    from scipy.ndimage import binary_erosion
    from scipy.ndimage.filters import gaussian_filter
    import numpy as np
    import os.path as op
    import math

    msk = nib.load( in_file )
    mskdata = msk.get_data()
    imshape = msk.get_shape()
    boundary = binary_erosion( mskdata ).astype(np.dtype('u1'))
    boundary = (( mskdata - boundary )).astype(float)
    
    x = ( int(0.4*imshape[0]), int(0.6*imshape[0]) )
    y = ( int(0.4*imshape[1]), int( imshape[1] ) )

    slice1 = np.zeros( shape=imshape )
    slice1[x[0]:x[1],y[0]:y[1],int(0.65*imshape[2])]  = 0.8

    slice2 = np.zeros( shape=imshape )
    slice2[x[0]:x[1],y[0]:y[1],int(0.35*imshape[2])]  = -0.6
    distortionfront = boundary * ( slice1 + slice2 )

    phasefield = gaussian_filter( distortionfront, sigma=5 )

    maxval = np.amax(phasefield)
    minval = np.amin(phasefield)

    phasefield = np.pi * (((phasefield + math.fabs(minval) ) / (maxval-minval)*2) - 1)
    phasefield = phasefield * 0.0005
    
    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_ph_unwrap_rads.nii.gz' % fname)

    nib.save( nib.Nifti1Image( phasefield, msk.get_affine(), msk.get_header() ), out_file )
    
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
