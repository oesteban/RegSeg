#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban - code@oscaresteban.es
# @Date:   2014-04-04 19:39:38
# @Last Modified by:   oesteban
# @Last Modified time: 2014-10-15 00:19:30

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from shutil import copyfileobj
import os
import os.path as op
import glob
import sys

import pyacwereg.workflows.evaluation as ev
from pyacwereg.workflows.smri import prepare_smri
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu

if __name__ == '__main__':
    parser = ArgumentParser(description='Run evaluation workflow',
                            formatter_class=RawTextHelpFormatter)

    g_input = parser.add_argument_group('Inputs')

    g_input.add_argument('-D', '--data_dir', action='store',
                         default=os.getenv('IXI_DATASET_HOME', os.getcwd()),
                         help='directory where subjects are found')
    g_input.add_argument('-s', '--subjects', action='store',
                         default='S*', help='subject id or pattern of ids')
    g_input.add_argument('-g', '--grid_size', action='store',
                         default=[6, 6, 6], nargs='+',
                         help='number of control points')
    g_input.add_argument('-w', '--work_dir', action='store',
                         default=os.getcwd(),
                         help='directory where subjects are found')
    g_input.add_argument('-N', '--name', action='store', default='EvaluationTests',
                         help='default workflow name, it will create a new folder')
    g_output = parser.add_argument_group('Outputs')
    g_output.add_argument('-o', '--out_csv', action='store',
                          help='output summary csv file')

    options = parser.parse_args()

    if not op.exists(options.work_dir):
        os.makedirs(options.work_dir)

    subjects_dir = op.join(options.data_dir, 'subjects')
    freesurfer_dir = op.join(options.data_dir, 'FREESURFER')

    sub_list = glob.glob(op.join(freesurfer_dir, options.subjects))
    subjects = [op.basename(sub) for sub in sub_list]

    if not len(subjects):
        print 'No subject was found in %s' % options.data_dir
        sys.exit(1)

    wf = pe.Workflow(name=options.name)
    wf.base_dir = options.work_dir
    infosource = pe.Node(niu.IdentityInterface(fields=['subject_id', 'data_dir']),
                         name="infosource")
    infosource.iterables = [('subject_id', subjects[0:3])]
    infosource.inputs.data_dir = options.data_dir

    prep = prepare_smri()
    bs = ev.bspline()
    bs.inputs.inputnode.grid_size = options.grid_size

    wf.connect([
        (infosource,  prep, [('subject_id', 'inputnode.subject_id'),
                             ('data_dir', 'inputnode.data_dir')]),
        (infosource,    bs, [('subject_id', 'inputnode.subject_id')]),
        (prep,          bs, [('outputnode.out_smri_brain', 'inputnode.in_file'),
                             ('outputnode.out_surfs', 'inputnode.in_surfs'),
                             ('outputnode.out_tpms', 'inputnode.in_tpms'),
                             ('outputnode.out_mask', 'inputnode.in_mask')])
    ])

    if options.out_csv is None:
        bs.inputs.inputnode.out_csv = op.join(
            options.work_dir, wf.name, 'results.csv')
    else:
        bs.inputs.inputnode.out_csv = options.out_csv

    try:
        wf.write_graph(
            graph2use='hierarchical', format='pdf', simple_form=True)
    except RuntimeError as e:
        print e

    try:
        wf.run()
    except RuntimeError as e:
        print e
