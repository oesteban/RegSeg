# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import os
import os.path as op

from nipype.interfaces.base import (traits, TraitedSpec, File, Directory,
                                    Undefined, isdefined, OutputMultiPath,
    InputMultiPath, BaseInterface, BaseInterfaceInputSpec)
from nipype.interfaces.io import IOBase


class FSFilesInputSpect( BaseInterfaceInputSpec ):
    fs_subjects_dir = Directory( os.environ['SUBJECTS_DIR'], exists=True, usedefault=True,
                                 desc='Freesurfer subjects dir' )
    subject_id = traits.String( mandatory=True, desc='Subject in Freesurfer system')

class FSFilesOutputSpect( TraitedSpec ):
    aseg = File(exists=True, desc='Path to aseg.mgz file' )
    rawavg = File(exists=True, desc='Path to rawavg.mgz file' )
    orig = File(exists=True, desc='Path to orig.mgz file' )
    norm = File(exists=True, desc='Path to norm.mgz file' )
    brain = File(exists=True, desc='Path to brain.mgz file' )

class FSFiles( IOBase ):
    """ Dumb datasource to get Freesurfer's files.

    Returns the path of different files of interest within the Freesurfer folder of a
    certain subject.

    Examples
    --------

    >>> from pyacwereg.utils.freesurfer import FSFiles
    >>> fs_files = FSFiles()
    >>> fs_files.inputs.fs_subjects_dir = './FREESURFER'
    >>> fs_files.inputs.subject_id = 'S400'
    >>> fs_files.run()
    """
    input_spec = FSFilesInputSpect
    output_spec = FSFilesOutputSpect

    def _list_outputs( self ):
        fs_dir = op.join( self.inputs.fs_subjects_dir, self.inputs.subject_id, 'mri' )
        outputs = self._outputs().get()
        outputs['aseg'] = op.join( fs_dir, 'aseg.mgz' )
        outputs['rawavg'] = op.join( fs_dir, 'rawavg.mgz' )
        outputs['orig'] = op.join( fs_dir, 'orig.mgz' )
        outputs['norm'] = op.join( fs_dir, 'norm.mgz' )
        outputs['brain'] = op.join( fs_dir, 'brain.mgz' )

        return outputs


def gen_fs_transform( in_file, fs_subjects_dir, subject_id, out_file=None ):
    """ Get transform matrix between orig (*tkRAS*) and native (*scannerRAS*)
    coordinates in Freesurfer. Implements first step of
    `this guide <http://surfer.nmr.mgh.harvard.edu/fswiki/FsAnat-to-NativeAnat>`_.

    Keyword arguments:
    in_file -- image in scanner coordinates
    fs_subjects_dir -- Freesurfer's SUBJECTS_DIR environment variable
    subject_id -- Subject identifier in Freesufer's folder tree
    out_file (optional) -- Filename for the transform matrix
    """
    import os.path as op
    import subprocess as sp

    in_orig = op.join( fs_subjects_dir, subject_id, 'mri', 'orig.mgz' )

    if out_file is None:
        out_file = op.abspath('./orig_to_native.dat' )

    sp.check_call( ["tkregister2", "--mov %s" % in_file, "--targ %s" % in_orig, "--reg %s" % out_file, "--noedit", "--regheader"])
    return out_file

def transform_surface( in_surf, in_reg, in_target, fs_subjects_dir, subject_id, out_file=None ):
    """ Apply a transform matrix to in_surf

    Keyword arguments:
    in_surf -- surface to be *resampled*
    in_reg -- Freesurfer's transform file
    in_target -- target image to get coordinates from
    fs_subjects_dir -- Freesurfer's SUBJECTS_DIR environment variable
    subject_id -- Subject identifier in Freesufer's folder tree
    out_file (optional) -- Filename for the resampled surface
    """
    import os
    import os.path as op
    import subprocess as sp

    my_env = os.environ.copy()
    my_env["SUBJECTS_DIR"] = "'%s'" % fs_subjects_dir

    hemi, surf = op.splitext(op.basename(in_surf))
    surf = surf[1:]

    if out_file is None:
        out_file = op.abspath('./%s.%s.native' % (hemi,surf) )

    try:
        cmd = "SUBJECTS_DIR='%s'; mri_surf2surf --sval-xyz %s --reg %s %s --tval %s --tval-xyz --hemi %s --s %s" % ( fs_subjects_dir, surf, in_reg, in_target, out_file, hemi, subject_id)
        sp.check_call( cmd, shell=True )
    except sp.CalledProcessError as e:
        print cmd
        raise

    return out_file

def MRIPreTess( in_file, in_norm, label=1, out_file="" ):
    """ Wraps mri_pretess command to generate a surface """
    import subprocess as sp
    import os
    import os.path as op

    if out_file=="":
        fname,ext = op.splitext( op.basename(in_file) )
        if ext==".gz":
            fname,_ = op.splitext(fname)

        out_file= op.join( op.abspath( os.getcwd() ), "%s_pretess.mgz" % fname )

    cmd="mri_pretess %s %d %s %s" % ( in_file, label, in_norm, out_file )
    proc = sp.check_call( cmd, shell=True )

    return out_file


def merge_labels( in_file, labels ):
    """ Binarizes a segmentation (`in_file`) by picking up labels only """
    import nibabel as nb
    import numpy as np
    import os.path as op
    from scipy import ndimage

    img = nb.load( in_file )
    imgdata = img.get_data()

    outdata = np.zeros( shape=imgdata.shape )
    for l in labels:
        outdata[ imgdata==l ] = 1

    #disk = disk_structure( (2,2,2) )
    #outdata = ndimage.binary_opening( outdata, structure=disk ).astype(np.int)
    #outdata = ndimage.binary_closing( outdata, structure=disk ).astype(np.int)

    fname, ext = op.splitext( in_file )

    if ext=='.gz':
        fname,_ = op.splitext(fname)

    out_file = '%s_bin.nii.gz' % fname
    nb.save( nb.Nifti1Image( outdata, img.get_affine(), img.get_header() ), out_file)
    return out_file

def fixvtk( in_file, in_ref, out_file=None ):
    """ Transforms an vtk file from Freesurfer's *tkRAS* coordinates
    to a target image in *scannerRAS* coordinates.
    """
    import nibabel as nb
    import numpy as np
    import os.path as op
    import subprocess as sp

    if out_file is None:
        fname,ext = op.splitext( op.basename(in_file) )
        if ext==".gz":
            fname,_ = op.splitext(fname)

        out_file= op.abspath( "%s_fixed.vtk" % fname )

    ref = nb.load( in_ref )
    cmd_info = "mri_info --tkr2scanner %s" % in_ref
    proc = sp.Popen( cmd_info, stdout=sp.PIPE, shell=True )
    data = bytearray( proc.stdout.read() )
    mstring =  np.fromstring(data.decode("utf-8"), sep='\n' )
    matrix = np.reshape(mstring, (4,-1) )

    with open( in_file, 'r' ) as f:
        with open( out_file, 'w+' ) as w:
            npoints = 0
            pointid = -5

            for i,l in enumerate(f):
                if (i==4):
                    s = l.split()
                    npoints = int( s[1] )
                    fmt = np.dtype( s[2] )
                elif( i>4 and pointid<npoints):
                    vert = [float(x) for x in l.split()]
                    vert.append( 1.0 )
                    newvert = np.dot( matrix, vert )
                    l = '%.9f  %.9f  %.9f\n' % tuple(newvert[0:3])

                w.write( l )
                pointid = pointid + 1

    return out_file
