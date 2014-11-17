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
            ('T1w/T2w_acpc_dc_restore.nii.gz', 'T2w_acpc_dc_restore.nii.gz')]

if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    from tmpfile import mkdtemp
    import subprocess as sp
    import shutil as sh

    curdir = os.getcwd()
    parser = ArgumentParser(description='Extract HCP data for TMI2015 paper',
                            formatter_class=RawTextHelpFormatter)
    g_input = parser.add_argument_group('Input')
    g_input.add_argument('-i', '--input', action='store')
    g_input.add_argument('-d', '--dest', action='store', default='.')

    opts = parser.parse_args()
    if opts.dest == '.':
        dest = os.getcwd()
    else:
        dest = op.abspath(opts.dest)

    for iname in opts.input:
        if not op.exists(iname):
            raise IOError('Input package not found')

        bname = op.basename(iname)
        caseid = bname.split('_', 1)[0]

        stname = op.dirname(iname) + bname.split('Diffusion', 1)[0]
        stname += 'Structural_preproc.zip'

        tmpdir = mkdtemp()
        os.chdir(tmpdir)
        try:
            sp.check_output(['unzip', op.abspath(iname)])
            sp.check_output(['unzip', op.abspath(stname)])

            dstfld = op.join(dest, 'HCP%s' % caseid)
            os.mkdir(dstfld)

            for fo, fd in filesmap:
                sh.copyfile(op.abspath(op.join(tmpdir, fo)),
                            op.abspath(dstfld, fd))
        except:
            print 'Extracting %s failed' % iname
            pass

        os.chdir(curdir)
        sh.rmtree(tmpdir)
