#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban - code@oscaresteban.es
# @Date:   2014-04-04 19:39:38
# @Last Modified by:   oesteban
# @Last Modified time: 2015-01-15 15:34:39


import os

if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    import os.path as op
    from glob import glob
    import numpy as np

    try:
        from enthought.etsconfig.api import ETSConfig
        ETSConfig.toolkit = 'null'
    except:
        pass

    from pyacwereg.workflows.realdata import hcp_workflow

    parser = ArgumentParser(description='PyACWEReg - Experiment on HCP data',
                            formatter_class=RawTextHelpFormatter)

    g_input = parser.add_argument_group('Input')
    g_input.add_argument('-S', '--subjects_dir', action='store',
                         default=os.getenv('NEURO_DATA_HOME', '..'),
                         help='directory where subjects should be found')
    g_input.add_argument('-s', '--subject', action='store', default='*',
                         nargs='+', help='subject id or pattern')
    g_input.add_argument('-w', '--work_dir', action='store',
                         default=os.getcwd(),
                         help='directory to store intermediate results')
    g_input.add_argument('-N', '--name', action='store',
                         default='TMI2015_EXP2',
                         help=('default workflow name, '
                               'it will create a new folder'))

    g_output = parser.add_argument_group('Output')
    g_output.add_argument('--out_csv', action='store',
                          help=('default output csv file'))

    g_options = parser.add_argument_group('Settings')
    g_options.add_argument('-T', '--type', action='store',
                           choices=['fmb', 'peb', 'fsl'], default='fmb',
                           help='select SDC workflow type')

    opts = parser.parse_args()

    if not op.exists(opts.work_dir):
        os.makedirs(opts.work_dir)

    data_dir = op.abspath(opts.subjects_dir)

    settings = {}
    settings['work_dir'] = opts.work_dir
    settings['data_dir'] = data_dir

    subj_list = []

    for subj in np.atleast_1d(opts.subject).tolist():
        subj = op.basename(subj)

        if '*' in subj:
            for subxpanded in glob(subj):
                subj_list.append(op.basename(subxpanded))
        else:
            subj_list.append(subj)

    settings['subject_id'] = subj_list

    if len(subj_list) == 0:
        raise RuntimeError('No subjects found in list')

    if opts.out_csv is None:
        settings['out_csv'] = op.join(opts.work_dir, opts.name, 'results.csv')
    else:
        settings['out_csv'] = opts.out_csv

    wf = hcp_workflow(name=opts.name, settings=settings)
    wf.base_dir = settings['work_dir']
    wf.write_graph(format='pdf')
    wf.run()
