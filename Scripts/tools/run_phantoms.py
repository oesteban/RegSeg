#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban - code@oscaresteban.es
# @Date:   2014-04-15 10:09:24
# @Last Modified by:   oesteban
# @Last Modified time: 2014-10-27 18:51:37

__author__ = "Oscar Esteban"
__copyright__ = "Copyright 2013, Biomedical Image Technologies (BIT), \
                 Universidad Politécnica de Madrid"
__credits__ = "Oscar Esteban"
__license__ = "FreeBSD"
__version__ = "0.1"
__maintainer__ = "Oscar Esteban"
__email__ = "code@oscaresteban.es"
__status__ = "Prototype"

try:
    from enthought.etsconfig.api import ETSConfig
    ETSConfig.toolkit = 'null'
except:
    pass
import os
import os.path as op


def phantoms_wf(options):
    import glob
    import nipype.pipeline.engine as pe
    from nipype.interfaces import utility as niu
    from pyacwereg.workflows import evaluation as ev

    subject_id = options.subject_id
    subject_dir = op.join(options.data_dir, subject_id)

    bs = ev.bspline(name=options.name)

    grid_size = options.grid_size
    if len(grid_size) == 1:
        grid_size = grid_size * 3

    bs.inputs.inputnode.grid_size = grid_size
    bs.inputs.inputnode.subject_id = subject_id

    if options.out_csv is None:
        bs.inputs.inputnode.out_csv = op.join(
            options.work_dir, bs.name, 'results.csv')
    else:
        bs.inputs.inputnode.out_csv = options.out_csv

    return bs

if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    from shutil import copyfileobj

    parser = ArgumentParser(description='Run evaluation workflow',
                            formatter_class=RawTextHelpFormatter)

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
    wf.base_dir = options.work_dir
    wf.write_graph(graph2use='hierarchical', format='pdf', simple_form=True)
    wf.run()
