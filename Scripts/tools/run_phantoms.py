#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban - code@oscaresteban.es
# @Date:   2014-04-15 10:09:24
# @Last Modified by:   oesteban
# @Last Modified time: 2014-10-22 10:09:27
import os.path as op


def phantoms_wf(options):
    import glob
    import nipype.pipeline.engine as pe
    from nipype.interfaces import utility as niu
    from pyacwereg.workflows import evaluation as ev

    subject_id = options.subject_id
    subject_dir = op.join(options.data_dir, subject_id)
    tpms = glob.glob(op.join(subject_dir, 'tpms', '*.nii.gz'))

    wf = pe.Workflow(name=options.name)
    wf.base_dir = options.work_dir
    infosource = pe.Node(niu.IdentityInterface(
        fields=['subject_id', 'data_dir', 'in_file', 'in_surfs', 'in_tpms',
                'in_mask']),
        name="infosource")

    infosource.inputs.subject_id = subject_id
    infosource.inputs.in_file = [
        op.join(subject_dir, 't1_weighted_SNR%d.nii.gz' % options.noise_snr),
        op.join(subject_dir, 't2_weighted_SNR%d.nii.gz' % options.noise_snr)]
    infosource.inputs.in_tpms = tpms
    infosource.inputs.in_surfs = glob.glob(
        op.join(subject_dir, 'surfs', '*.vtk'))
    infosource.inputs.in_mask = op.join(subject_dir, 'signal_mask.nii.gz')

    bs = ev.bspline(n_tissues=len(tpms))
    bs.inputs.inputnode.grid_size = options.grid_size

    wf.connect([
        (infosource,    bs, [('subject_id', 'inputnode.subject_id'),
                             ('in_file', 'inputnode.in_file'),
                             ('in_surfs', 'inputnode.in_surfs'),
                             ('in_tpms', 'inputnode.in_tpms'),
                             ('in_mask', 'inputnode.in_mask')])
    ])

    if options.out_csv is None:
        bs.inputs.inputnode.out_csv = op.join(
            options.work_dir, wf.name, 'results.csv')
    else:
        bs.inputs.inputnode.out_csv = options.out_csv

    return wf

if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    from shutil import copyfileobj
    import os
    parser = ArgumentParser(description='Run evaluation workflow',
                            formatter_class=RawTextHelpFormatter)
    try:
        from enthought.etsconfig.api import ETSConfig
        ETSConfig.toolkit = 'null'
    except:
        pass

    g_input = parser.add_argument_group('Inputs')

    g_input.add_argument(
        '-D', '--data_dir', action='store',
        default=op.join(os.getenv('NEURO_DATA_HOME', os.getcwd()), 'phantoms'),
        help='directory where subjects are found')
    g_input.add_argument(
        '-s', '--subject_id', action='store', default='S001',
        help='selects phantom\'s shape model')
    g_input.add_argument(
        '-n', '--noise_snr', action='store', default=30, type=int,
        help='selects phantom\'s number of regions')
    g_input.add_argument(
        '-g', '--grid_size', action='store', default=[6, 6, 6], nargs='+',
        type=int, help='number of control points')
    g_input.add_argument(
        '-w', '--work_dir', action='store', default=os.getcwd(),
        help='directory where subjects are found')
    g_input.add_argument(
        '-N', '--name', action='store', default='PhantomTests',
        help='default workflow name, it will create a new folder')
    g_output = parser.add_argument_group('Outputs')
    g_output.add_argument(
        '-o', '--out_csv', action='store', help='output summary csv file')

    options = parser.parse_args()

    if not op.exists(options.work_dir):
        os.makedirs(options.work_dir)

    wf = phantoms_wf(options)
    wf.write_graph(graph2use='hierarchical', format='pdf', simple_form=True)
    wf.run()
