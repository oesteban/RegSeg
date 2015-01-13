#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban - code@oscaresteban.es
# @Date:   2014-04-15 10:09:24
# @Last Modified by:   oesteban
# @Last Modified time: 2015-01-13 12:21:46

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

    subject_id = '%s_snr%03d' % (options.shape, options.snr)
    subject_dir = op.join(options.data_dir, subject_id)

    bs = ev.bspline(name=options.name)

    grid_size = options.grid_size
    if len(grid_size) == 1:
        grid_size = grid_size * 3

    bs.inputs.inputnode.grid_size = grid_size
    bs.inputs.inputnode.subject_id = subject_id
    bs.inputs.inputnode.lo_matrix = options.lo_matrix
    bs.inputs.inputnode.hi_matrix = options.hi_matrix
    bs.inputs.inputnode.shape = options.shape
    bs.inputs.inputnode.snr = options.snr
    bs.inputs.inputnode.cortex = options.no_cortex

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
        '-s', '--shape', action='store', default='gyrus',
        help='selects phantom\'s shape model')
    g_input.add_argument(
        '-n', '--snr', action='store', default=400, type=int,
        help='generate signal with certain SNR')
    g_input.add_argument(
        '--no_cortex', action='store_false',
        help='do not generate cortex-like crust')

    g_input.add_argument(
        '--lo_matrix', action='store', default=51, type=int,
        help='low-resolution matrix size')
    g_input.add_argument(
        '--hi_matrix', action='store', default=101, type=int,
        help='hi-resolution matrix size')

    g_input.add_argument(
        '-g', '--grid_size', action='store', default=[4, 4, 4], nargs='+',
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
