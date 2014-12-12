#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2014-12-11 15:08:23
# @Last Modified by:   oesteban
# @Last Modified time: 2014-12-12 19:20:25


def plot_report(df, out_file=None):
    import nibabel as nb
    import numpy as np
    import pandas as pd
    import seaborn as sn
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as mgs

    plt.style.use('ggplot')
    c = plt.rcParams['axes.color_cycle'][0]
    levels = pd.Series(df.level).unique()

    fd = {'size': 20, 'color': 'w', 'weight': 'bold'}
    lposx = .05
    lposy = .9
    lha = 'left'

    fig = plt.figure(figsize=(18, 19))
    outergs = mgs.GridSpec(2, 2, wspace=0.2, hspace=0.2)

    innergs0 = mgs.GridSpecFromSubplotSpec(
        2, len(levels), subplot_spec=outergs[0], hspace=0.3, wspace=0.08)
    allax1 = []
    allax2 = []
    for i, l in enumerate(levels):
        ax1 = plt.Subplot(fig, innergs0[i])
        allax1.append(ax1)
        ax2 = plt.Subplot(fig, innergs0[i + len(levels)])
        allax2.append(ax2)
        ldf = df[df.level == i]
        ldf.plot(ax=ax1, x='iteration', y='E_t')
        fig.add_subplot(ax1)
        ldf.plot(ax=ax2, x='iteration', y='step')
        fig.add_subplot(ax2)
        ax1.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax1.transAxes, fontdict=fd)
        ax2.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax2.transAxes, fontdict=fd)

    innergs3 = mgs.GridSpecFromSubplotSpec(
        2, len(levels), subplot_spec=outergs[2], hspace=0.3, wspace=0.05)
    allax1 = []
    allax2 = []
    for i, l in enumerate(levels):
        ldf = df[df.level == i]
        ax1 = plt.Subplot(fig, innergs3[i])
        allax1.append(ax1)
        ax2 = plt.Subplot(fig, innergs3[i + len(levels)])
        allax2.append(ax2)

        ldf.plot(ax=ax1, x='iteration', y=['norm', 'max_gk'],
                 secondary_y=['norm'])
        # ldf.plot(ax=ax1, x='iteration', y='max_gk', secondary_y=True)
        fig.add_subplot(ax1)

        ldf.plot(ax=ax2, x='iteration', y='g_50')
        ldf.plot(ax=ax2, x='iteration', y='g_05', color=c, style='--', lw=.7)
        ldf.plot(ax=ax2, x='iteration', y='g_95', color=c, style='--', lw=.7)
        ax2.fill_between(
            ldf.iteration, ldf['g_25'], ldf['g_75'], alpha=0.2, color=c)
        ax2.fill_between(
            ldf.iteration, ldf['g_05'], ldf['g_95'], alpha=0.05, color=c)
        fig.add_subplot(ax2)
        ax1.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax1.transAxes, fontdict=fd)
        ax2.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax2.transAxes, fontdict=fd)

    innergs1 = mgs.GridSpecFromSubplotSpec(
        1, len(levels), subplot_spec=outergs[1], hspace=0.0, wspace=0.05)
    for i, l in enumerate(levels):
        ldf = df[df.level == i]
        ax1 = plt.Subplot(fig, innergs1[i])

        eid = 0
        erois = []
        while 'E_%02d' % eid in ldf:
            erois.append('E_%02d' % eid)
            eid += 1

        ldf.plot(ax=ax1, x='iteration', y=erois, kind='area', stacked=True)
        fig.add_subplot(ax1)
        ax1.text(.9, .05, 'Level %d' %
                 (i + 1), ha='right', transform=ax1.transAxes, fontdict=fd)

    legend = ax1.legend(loc='center right')
    frame = legend.get_frame()
    frame.set_facecolor('w')

    innergs2 = mgs.GridSpecFromSubplotSpec(
        1, len(levels), subplot_spec=outergs[3], hspace=0.0, wspace=0.05)
    for i, l in enumerate(levels):
        ax = plt.Subplot(fig, innergs2[i])
        ldf.plot(ax=ax, x='iteration', y=['E_d', 'E_r'],
                 kind='area', stacked=True)
        plt.plot(ldf['E_t'])
        fig.add_subplot(ax)
        ax.text(.9, .05, 'Level %d' %
                (i + 1), ha='right', transform=ax.transAxes, fontdict=fd)

    legend = ax.legend(loc='upper right', bbox_to_anchor=(0.0, -0.05), ncol=3)
    frame = legend.get_frame()
    frame.set_facecolor('w')

    fig.suptitle('Convergence report', fontsize=40)

    gps = outergs.get_grid_positions(fig)
    fig.text(gps[3][0] - 0.38, gps[3][1] + 0.02,
             'Total energy evolution', size=20)
    fig.text(gps[3][0] - 0.38, gps[3][1] - 0.185,
             r'Step size ($\delta$) evolution', size=20)
    fig.text(gps[0][0], gps[3][1] + 0.02,
             'Region-wise evolution of energy', size=20)
    fig.text(gps[3][0] - 0.38, gps[1][1] + 0.02,
             r'Update: $\max_{k}\{\Vert \mathbf{u}_k \Vert\}$', size=20)
    fig.text(gps[3][0] - 0.38, gps[1][1] - 0.185,
             r'Gradient at vertices ($\mathbf{g}_i$) distribution', size=20)
    fig.text(gps[0][0], gps[1][1] + 0.02, 'Energy decomposition', size=20)

    if out_file is not None:
        plt.savefig(out_file, format='pdf', dpi=300, bbox_inches='tight')
    else:
        plt.show()
