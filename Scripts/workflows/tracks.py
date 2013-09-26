#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os, os.path as op
import warnings
import numpy as np
import nibabel as nib
from nipype.interfaces.base import (TraitedSpec, File, InputMultiPath,
                                    OutputMultiPath, Undefined, traits,
                                    isdefined, OutputMultiPath, 
                                    CommandLineInputSpec, CommandLine,
                                    BaseInterface, BaseInterfaceInputSpec,
                                    traits )
from nipype.utils.filemanip import split_filename,fname_presuffix



class ConnectivityMatrixInputSpec(BaseInterfaceInputSpec):
    in_file = File( exists=True, mandatory=True,
                    desc='streamline TrackVis trk file' )
    roi_file= File( exists=True, mandatory=True,
                    desc='file containing ROIs to create nodes')
    min_length = traits.Float( 10.0, usedefault=True,
                               desc='minimum length (mm.) of tracks to be considered' )
    out_prefix = File( desc='path prefix that will be used for output files' )
    
class ConnectivityMatrixOutputSpec(TraitedSpec):
    out_conmat = File( desc='textfile containing connectivity tracks count' )
    out_lenmat = File( desc='textfile containing the average length of tracks per connection' )
    out_tracks = File( desc='TrackVis trk file with streamlines considered in output matrix' )

class ConnectivityMatrix(BaseInterface):
    """
    Builds a connectivity matrix from a file containing seeding ROIs and a 
    trackvis format streamline file.

    Example
    -------

    >>> import tracks as trks
    >>> intf = trks.ConnectivityMatrix()
    >>> intf.inputs.in_file = trk_filename
    >>> intf.inputs.roi_file = roi_filename
    >>> intf.inputs.out_prefix = op.join( working_dir, 'interfacetest' )
    >>> res = intf.run()
    """

    input_spec = ConnectivityMatrixInputSpec
    output_spec = ConnectivityMatrixOutputSpec

    def _run_interface(self, runtime):
        trkfile, hdr = nib.trackvis.read( self.inputs.in_file )
        roifile = nib.load( self.inputs.roi_file )
        roidata = roifile.get_data()
        nROI = np.max(roidata.reshape(-1))
        min_length = self.inputs.min_length

        out_trk, out_ft = filter_track( trkfile, roidata, trkhdr=hdr, min_length=min_length )
        
        conmat = np.zeros((nROI,nROI), dtype='uint16')
        lenmat = np.zeros((nROI,nROI), dtype='f32')
    
        for t in out_ft:
            start = t[0][0]-1
            end   = t[0][1]-1
            l = t[1]
            
            if start>end:
                (start,end) = (end,start)
                
            conmat[start,end] = conmat[start,end] + 1
            lenmat[start,end] = lenmat[start,end] + l
        
        lenmat[lenmat>0] = lenmat[lenmat>0]/conmat.astype( 'f32' )[lenmat>0]
        
        conmat = conmat + conmat.T
        lenmat = lenmat + lenmat.T

        if not isdefined( self.inputs.out_prefix ):
            fname,fext = op.splitext( op.basename( self.inputs.in_file ) )
            if fext=='.gz':
                fname,fext2 = op.splitext(fname)
                fext = fext2+fext
            
            self.inputs.out_prefix = op.abspath( fname )
        
        nib.trackvis.write( self.inputs.out_prefix + '_tracks.trk', out_trk, hdr_mapping=hdr )
        np.savetxt( self.inputs.out_prefix + '_conmat.txt', conmat, fmt='%d' )
        np.savetxt( self.inputs.out_prefix + '_lenmat.txt', lenmat, fmt='%.3f' )
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["out_conmat"] = self.inputs.out_prefix + '_conmat.txt'
        outputs["out_lenmat"] = self.inputs.out_prefix + '_lenmat.txt'
        outputs["out_tracks"] = self.inputs.out_prefix + '_tracks.trk'
        return outputs


def track_analyze( points, roidata, vox_size=(1.0,1.0,1.0) ):
    startval = 0
    startidx = -1
    vox_size=np.array( vox_size, dtype='f32' )
    im_idx = np.array( np.shape(roidata) )-1

    invalid_idx = []
    for i,point in enumerate(points):
        vox = tuple( ( np.around( point / vox_size ) ).astype( 'uint16' ))
        if np.any( vox > im_idx ):
            invalid_idx+= [i]    

    points = np.delete( points, invalid_idx, 0 )    
    nPoints = np.shape(points)[0]

    if nPoints == 0:
        return (None,0)

    for i,point in enumerate(points):
        vox = tuple( ( np.around( point / vox_size ) ).astype( 'uint16' ))
        val = roidata[vox]
        if val>0:
            startval=val
            startidx = i
            break
            
            
    endval = 0
    endidx = -1
    for i,point in enumerate(reversed(points)):
        vox = tuple( ( np.around( point / vox_size ) ).astype( 'uint16' ))
        val = roidata[vox]
        if val>0:
            endval=val
            endidx = nPoints-1-i
            break
    
    length = 0.0
    
    new_track = []
    
    for i,point in enumerate(points):
        if i>= startidx:
            new_track.append(point)
            if i<endidx:
                length+= np.linalg.norm( point - points[i+1] )
            if i==endidx:
                break
                   
    track_features = [ (startval,endval), length ]
    
    return (new_track, track_features)

def filter_track( trkfile, roidata, trkhdr=None, min_length=3.0 ):
    streams = [ ss[0] for ss in trkfile ]
    
    if not trkhdr is None:
        vox_size = trkhdr['voxel_size']
    
    filtered_ft  = []
    out_trk = []
    
    for idx,trk in enumerate(streams):
        nt,tf = track_analyze( trk, roidata, vox_size )

        if nt is None:
            continue
       
        if tf[-1]>min_length:
            if tf[0][0]!=tf[0][1]:
                filtered_ft.append( tf )
                new_trk = (np.array(nt), trkfile[idx][1],trkfile[idx][2])
                out_trk.append(new_trk)
            
    return out_trk, filtered_ft


