#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2014-11-17 15:32:49
# @Last Modified by:   oesteban
# @Last Modified time: 2014-11-17 16:13:07

import os
import os.path as op

filesmap = [('T1w/aparc+aseg.nii.gz', 'aparc+aseg.nii.gz'),
            ('T1w/brainmask_fs.nii.gz', 'brainmask_fs.nii.gz'),
            ('T1w/Diffusion/bvals', 'bvals'),
            ('T1w/Diffusion/bvecs', 'bvecs'),
            ('T1w/Diffusion/data.nii.gz', 'dwi.nii.gz'),
            ('T1w/Diffusion/nodif_brain_mask.nii.gz', 'dwi_brainmask.nii.gz'),
            ('T1w/T1w_acpc_dc_restore.nii.gz', 'T1w_acpc_dc_restore.nii.gz'),
            ('T1w/T2w_acpc_dc_restore.nii.gz', 'T2w_acpc_dc_restore.nii.gz'),
            ('T1w/T1w_acpc_dc_restore_brain.nii.gz',
             'T1w_acpc_dc_restore_brain.nii.gz'),
            ('T1w/T2w_acpc_dc_restore_brain.nii.gz',
             'T2w_acpc_dc_restore_brain.nii.gz')]

if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    from tempfile import mkdtemp
    import subprocess as sp
    import shutil as sh

    curdir = os.getcwd()
    parser = ArgumentParser(description='Extract HCP data for TMI2015 paper',
                            formatter_class=RawTextHelpFormatter)
    g_input = parser.add_argument_group('Input')
    g_input.add_argument('-i', '--input', nargs='+', action='store')
    g_input.add_argument('-d', '--dest', action='store', default='.')

    opts = parser.parse_args()
    if opts.dest == '.':
        dest = os.getcwd()
    else:
        dest = op.abspath(opts.dest)

    for iname in opts.input:
        if not op.exists(iname):
            raise IOError('Input package not found')

        iname = op.abspath(iname)

        bname = op.basename(iname)
        dname = op.dirname(iname)
        caseid = bname.split('_', 1)[0]

        stname = op.join(dname, bname.split('Diffusion', 1)[0])
        stname += 'Structural_preproc.zip'

        tmpdir = mkdtemp()
        os.chdir(tmpdir)
        try:
            sp.check_output(['unzip', iname])
            sp.check_output(['unzip', stname])
            dstfld = op.join(dest, 'HCP%s' % caseid)

            try:
                sh.rmtree(dstfld)
            except:
                pass

            os.mkdir(dstfld)
            
            for f in filesmap:
                forig = op.abspath(op.join(tmpdir, caseid, f[0]))
                fdest = op.abspath(op.join(dstfld, f[1]))
                print 'Copying: %s -> %s' % (forig, fdest)
                sh.copy(forig, fdest)
        except Exception as e:
            print 'Extracting %s failed' % iname
            print e
            print os.listdir(tmpdir)
            pass

        os.chdir(curdir)
        sh.rmtree(tmpdir)
