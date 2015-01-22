#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# @Author: oesteban - code@oscaresteban.es
# @Date:   2014-04-15 10:09:24
# @Last Modified by:   oesteban
# @Last Modified time: 2015-01-22 13:00:33

try:
    from enthought.etsconfig.api import ETSConfig
    ETSConfig.toolkit = 'null'
except:
    pass
import os
import os.path as op


def phantoms_wf(options, cfg):
    import glob
    import nipype.pipeline.engine as pe
    from nipype import config, logging
    from nipype.interfaces import utility as niu
    from pyacwereg.workflows import evaluation as ev

    config.update_config(cfg)
    logging.update_logging(config)

    grid_size = options.grid_size
    if len(grid_size) == 1:
        grid_size = grid_size * 3

    bs = ev.bspline(name=options.name, shapes=options.shape,
                    snr_list=options.snr,
                    N=options.repetitions)
    bs.inputs.inputnode.grid_size = grid_size
    bs.inputs.inputnode.lo_matrix = options.lo_matrix
    bs.inputs.inputnode.hi_matrix = options.hi_matrix
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

    # Experimental Settings
    ################################
    g_input = parser.add_argument_group('Experimental Settings')
    g_input.add_argument(
        '-s', '--shape', action='store', default='gyrus', nargs='+',
        help='selects phantom\'s shape model')
    g_input.add_argument(
        '-n', '--snr', action='store', default=400, nargs='+', type=int,
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

    g_input.add_argument('-R', '--repetitions', action='store', default=1,
                         type=int, help='number of repetitions')

    # General Settings
    ################################
    g_settings = parser.add_argument_group('General settings')
    g_settings.add_argument(
        '-w', '--work_dir', action='store', default=os.getcwd(),
        help='directory where subjects are found')
    g_settings.add_argument(
        '-N', '--name', action='store', default='PhantomTests',
        help='default workflow name, it will create a new folder')
    g_settings.add_argument('--nthreads', action='store', default=0,
                            type=int, help='number of repetitions')
    g_settings.add_argument(
        '-D', '--data_dir', action='store',
        default=op.join(os.getenv('NEURO_DATA_HOME', os.getcwd()), 'phantoms'),
        help='directory where subjects are found')
    g_settings.add_argument('--debug', action='store_true', default=False,
                            help='switch debug mode ON')

    # Outputs
    ################################
    g_output = parser.add_argument_group('Outputs')
    g_output.add_argument(
        '-o', '--out_csv', action='store', help='output summary csv file')

    options = parser.parse_args()

    # Setup multiprocessing
    nthreads = options.nthreads
    if nthreads == 0:
        from multiprocessing import cpu_count
        nthreads = cpu_count()

    cfg = {}
    cfg['plugin'] = 'Linear'
    if nthreads > 1:
        cfg['plugin'] = 'MultiProc'
        cfg['plugin_args'] = {'n_proc': nthreads}

    # Setup work_dir
    if not op.exists(options.work_dir):
        os.makedirs(options.work_dir)

    # Setup logging dir
    log_dir = op.abspath('logs')
    cfg['logging'] = {'log_directory': log_dir, 'log_to_file': True,
                      'workflow_level': 'INFO', 'interface_level': 'INFO'}
    if not op.exists(log_dir):
        os.makedirs(log_dir)

    # Setup debug mode
    if options.debug:
        cfg['logging']['workflow_level'] = 'DEBUG'
        cfg['logging']['interface_level'] = 'DEBUG'

    wf = phantoms_wf(options, cfg)
    wf.base_dir = options.work_dir
    wf.write_graph(graph2use='hierarchical', format='pdf', simple_form=True)
    wf.run()
