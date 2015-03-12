#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2014-11-17 15:32:49
# @Last Modified by:   oesteban
# @Last Modified time: 2015-02-10 12:13:50

filesmap = {}

filesmap['struct_preproc'] = [
    ('T1w/aparc+aseg.nii.gz', 'aparc+aseg.nii.gz'),
    ('T1w/brainmask_fs.nii.gz', 'brainmask_fs.nii.gz'),
    ('T1w/T1w_acpc_dc_restore.nii.gz', 'T1w_acpc_dc_restore.nii.gz'),
    ('T1w/T2w_acpc_dc_restore.nii.gz', 'T2w_acpc_dc_restore.nii.gz'),
    ('T1w/T1w_acpc_dc_restore_brain.nii.gz',
     'T1w_acpc_dc_restore_brain.nii.gz'),
    ('T1w/T2w_acpc_dc_restore_brain.nii.gz',
     'T2w_acpc_dc_restore_brain.nii.gz')]

filesmap['diff_preproc'] = [
    ('T1w/Diffusion/bvals', 'bvals'),
    ('T1w/Diffusion/bvecs', 'bvecs'),
    ('T1w/Diffusion/data.nii.gz', 'dwi.nii.gz'),
    ('T1w/Diffusion/nodif_brain_mask.nii.gz', 'dwi_brainmask.nii.gz')]


filesmap['struct_unproc'] = [
    ('unprocessed/{tesla}/T1w_MPR1/{caseid}_{tesla}_FieldMap_Magnitude.nii.gz',
     'FM_mag.nii.gz'),
    ('unprocessed/{tesla}/T1w_MPR1/{caseid}_{tesla}_FieldMap_Phase.nii.gz',
     'FM_pha.nii.gz')]

defparms = {'enc_dir': 'y-',
            'echospacing': 7.800117313764398e-4,
            'delta_te': 2.46e-3,
            'epi_factor': 128,
            'epi_lines': 57,
            'epi_acc': 2,
            'field_strength': 3.0,
            'field_axis': 'z'
            }

if __name__ == '__main__':
    import os
    import errno
    import os.path as op
    import json
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    from tempfile import mkdtemp
    import subprocess as sp
    import shutil as sh

    curdir = os.getcwd()
    parser = ArgumentParser(description='Extract HCP data for diffusion experiments',
                            formatter_class=RawTextHelpFormatter)
    g_input = parser.add_argument_group('Input')
    g_input.add_argument('-i', '--input', nargs='+', action='store')
    g_input.add_argument('-d', '--dest', action='store', default='.')
    g_input.add_argument('-t', '--tmpdir', action='store', default='/tmp/')

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
        bnameparts = bname.split('_')
        caseid = bnameparts[0]
        tesla = bnameparts[1]

        stname = op.join(dname, bname.split('Diffusion', 1)[0])
        stname += 'Structural_preproc.zip'

        stname2 = op.join(dname, bname.split('Diffusion', 1)[0])
        stname2 += 'Structural_unproc.zip'

        tmpdir = mkdtemp(prefix=opts.tmpdir)
        zipnames = (iname, stname, stname2)
        os.chdir(tmpdir)
        try:
            dstfld = op.join(dest, 'HCP%s' % caseid)
            try:
                os.makedirs(dstfld)
            except OSError as exc:  # Python >2.5
                if exc.errno == errno.EEXIST and op.isdir(dstfld):
                    pass
                else:
                    raise

            for n, v in zip(zipnames, ['diff_preproc',
                                       'struct_preproc', 'struct_unproc']):
                print 'Extracting HCP%s_%s to %s' % (caseid, tesla, tmpdir)
                ofiles = [f[0].format(tesla=tesla, caseid=caseid)
                          for f in filesmap[v]]
                dfiles = [f[1].format(tesla=tesla, caseid=caseid)
                          for f in filesmap[v]]

                rewrite = [fid for fid, f in enumerate(dfiles)
                           if op.exists(op.join(dstfld, f))]

                rfiles = []
                for idx in reversed(rewrite):
                    rfiles.append(ofiles.pop(idx))
                    dfiles.pop(idx)

                if len(ofiles) > 0:
                    flist = [op.join(caseid, f) for f in ofiles]
                    sp.check_output(['unzip', n] + flist)
                    for of, df in zip(flist, dfiles):
                        forig = op.abspath(of)
                        fdest = op.abspath(op.join(dstfld, df))
                        print 'Copying: %s -> %s' % (forig, fdest)
                        sh.copy(forig, fdest)

                if len(rfiles) > 0:
                    print 'Skipped %s' % ' '.join(rfiles)

            if not op.exists(op.join(dstfld, 'parameters.txt')):
                with open(op.join(dstfld, 'parameters.txt'), 'w') as pf:
                    json.dump(defparms, pf)

        except Exception as e:
            print 'Unpacking %s failed' % iname
            print e
            print os.listdir(tmpdir)
            pass

        os.chdir(curdir)
        sh.rmtree(tmpdir)
