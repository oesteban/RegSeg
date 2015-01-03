#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2014-12-11 15:08:23
# @Last Modified by:   oesteban
# @Last Modified time: 2015-01-03 13:25:32


def plot_report(df, levels_df=None, out_file=None):
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

    for i, l in enumerate(levels):
        ax1 = plt.Subplot(fig, innergs0[i])
        ax2 = plt.Subplot(fig, innergs0[i + len(levels)])
        fig.add_subplot(ax1)
        fig.add_subplot(ax2)

        ldf = df[df.level == i]
        plt.sca(ax1)
        ldf.E_t.plot(label=r'$E_T(t)$')

        its = ldf.iteration[ldf['desc_update'] != 0]
        if not its.empty:
            bottom = ax1.get_ylim()[0]
            dupdt = ldf.desc_update[its].astype(np.float32) * ldf.E_t[its]
            for it, v in zip(its, dupdt):
                plt.plot((it, it), (bottom, v), color=c, linewidth=0.7,
                         alpha=0.9)

        ldf.step.plot(secondary_y=True, style='k--', linewidth=0.7,
                      label=r'$\delta(t)$')
        ax1.set_xlabel('iteration')

        if levels_df is not None:
            if not levels_df.empty:
                th_e = levels_df['total'][i]
                ax1.plot([th_e] * (len(ldf['E_t'])))

        plt.sca(ax2)
        ldf.max_uk.plot(label=r'$\max_{k}\{\Vert \mathbf{u}_k(t) \Vert\}$')

        if not its.empty:
            bottom = ax2.get_ylim()[0]
            dupdt = ldf.desc_update[its].astype(np.float32) * ldf.max_uk[its]
            for it, v in zip(its, dupdt):
                plt.plot((it, it), (bottom, v), color=c, linewidth=0.7,
                         alpha=0.9)
        ldf.max_gk.plot(secondary_y=True,
                        label=r'$\max_{k}\{\Vert \mathbf{g}_k(t) \Vert\}$')

        ax2.set_xlabel('iteration')

        ax1.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax1.transAxes, fontdict=fd)
        ax2.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax2.transAxes, fontdict=fd)

    innergs3 = mgs.GridSpecFromSubplotSpec(
        1, len(levels), subplot_spec=outergs[2], hspace=0.3, wspace=0.05)

    for i, l in enumerate(levels):
        ax1 = plt.Subplot(fig, innergs3[i])
        fig.add_subplot(ax1)

        ldf = df[df.level == i]
        plt.sca(ax1)
        ldf.plot(x='iteration', y='g_50')
        ldf.plot(x='iteration', y='g_05', color=c, style='--', lw=.7)
        ldf.plot(x='iteration', y='g_95', color=c, style='--', lw=.7)
        ax1.fill_between(
            ldf.iteration, ldf['g_25'], ldf['g_75'], alpha=0.2, color=c)
        ax1.fill_between(
            ldf.iteration, ldf['g_05'], ldf['g_95'], alpha=0.05, color=c)
        ax1.set_xlabel('iteration')
        ax1.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax1.transAxes, fontdict=fd)

        # ax2 = plt.Subplot(fig, innergs3[i + len(levels)])
        # fig.add_subplot(ax2)
        # ldf.plot(ax=ax2, x='iteration', y='g_50')
        # ldf.plot(ax=ax2, x='iteration', y='g_05', color=c, style='--', lw=.7)
        # ldf.plot(ax=ax2, x='iteration', y='g_95', color=c, style='--', lw=.7)
        # ax2.fill_between(
        #     ldf.iteration, ldf['g_25'], ldf['g_75'], alpha=0.2, color=c)
        # ax2.fill_between(
        #     ldf.iteration, ldf['g_05'], ldf['g_95'], alpha=0.05, color=c)
        # ax2.text(lposx, lposy, 'Level %d' %
        #          (i + 1), ha=lha, transform=ax2.transAxes, fontdict=fd)

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
             'Update evolution', size=20)
    fig.text(gps[0][0], gps[3][1] + 0.02,
             'Region-wise evolution of energy', size=20)
    fig.text(gps[3][0] - 0.38, gps[1][1] + 0.02,
             r'Gradient at vertices ($\mathbf{g}_i$) distribution', size=20)

    fig.text(gps[0][0], gps[1][1] + 0.02, 'Energy decomposition', size=20)

    if out_file is not None:
        plt.savefig(out_file, format='pdf', dpi=300, bbox_inches='tight')
    else:
        plt.show()


def jointplot_data(im1data, im2data, in_seg, labels=None, out_file=None,
                   f1lims=None, f2lims=None,
                   f1name='Image A', f2name='Image B'):
    import os.path as op
    import nibabel as nb
    import numpy as np
    import os
    import seaborn as sn
    import pandas as pd
    import matplotlib.pyplot as plt

    if f1lims is None:
        f1lims = (im1data.min(), im1data.max())

    if f2lims is None:
        f2lims = (im2data.min(), im2data.max())

    df_pieces = []
    for i, l in enumerate(np.unique(in_seg)):
        filt = im1data[in_seg == l]
        d = {'tissue': [labels[i]] * len(filt)}
        d[f1name] = filt
        d[f2name] = im2data[in_seg == l]
        df_pieces.append(pd.DataFrame(d))
    df = pd.concat(df_pieces)
    sn.set_context("talk", font_scale=1.8)
    g = sn.JointGrid(f1name, f2name, df, hue='tissue',
                     xlim=f1lims, ylim=f2lims, size=20)
    g.plot_joint(sn.kdeplot, linewidths=3.0)
    g.plot_marginals(sn.distplot)

    nlabels = len(g.ax_joint.yaxis.get_ticklabels())

    plt.setp(g.ax_joint.get_yticklabels(), visible=False)
    plt.setp(g.ax_joint.get_xticklabels(), visible=False)

    if out_file is None:
        out_file = op.abspath('jointplot.pdf')

    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    return out_file


def jointplot_gmm(loc, cov, labels=None, out_file=None,
                  xlims=None, ylims=None,
                  xname='Image A', yname='Image B',
                  size=20, ratio=5, space=.2,):
    import os.path as op
    import nibabel as nb
    import numpy as np
    import os
    import seaborn as sn
    from seaborn.palettes import color_palette, light_palette
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
    from matplotlib import patches as mpatches

    if len(loc) != len(cov):
        raise RuntimeError('Mixture model does not contain elements')

    f = plt.figure(figsize=(size, size))
    gs = plt.GridSpec(ratio + 1, ratio + 1)
    sn.set_context("talk", font_scale=1.8)

    ax_joint = f.add_subplot(gs[1:, :-1])
    ax_marg_x = f.add_subplot(gs[0, :-1], sharex=ax_joint)
    ax_marg_y = f.add_subplot(gs[1:, -1], sharey=ax_joint)

    # Turn off tick visibility for the measure axis on the marginal plots
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    # Turn off the ticks on the density axis for the marginal plots
    plt.setp(ax_marg_x.yaxis.get_majorticklines(), visible=False)
    plt.setp(ax_marg_x.yaxis.get_minorticklines(), visible=False)
    plt.setp(ax_marg_y.xaxis.get_majorticklines(), visible=False)
    plt.setp(ax_marg_y.xaxis.get_minorticklines(), visible=False)
    plt.setp(ax_marg_x.get_yticklabels(), visible=False)
    plt.setp(ax_marg_y.get_xticklabels(), visible=False)
    plt.setp(ax_joint.get_yticklabels(), visible=False)
    plt.setp(ax_joint.get_xticklabels(), visible=False)

    ax_joint.set_xlabel(xname)
    ax_joint.set_ylabel(yname)

    ax_marg_x.yaxis.grid(False)
    ax_marg_y.xaxis.grid(False)
    palette = color_palette("husl", n_colors=len(loc))

    # Make the grid look nice
    sn.utils.despine(f)
    sn.utils.despine(ax=ax_marg_x, left=True)
    sn.utils.despine(ax=ax_marg_y, bottom=True)
    f.tight_layout()
    f.subplots_adjust(hspace=space, wspace=space)

    if xlims is None:
        xlims = (0.0, 1.0)

    if ylims is None:
        ylims = (0.0, 1.0)

    xstep = (xlims[1] - xlims[0]) / 100
    ystep = (ylims[1] - ylims[0]) / 100

    x = np.arange(xlims[0], xlims[1] + xstep, xstep)
    y = np.arange(ylims[0], ylims[1] + ystep, ystep)
    X, Y = np.meshgrid(x, y)

    if labels is None:
        labels = [None] * len(loc)

    patches = []
    for mu, sigma, color, l in zip(loc, cov, palette, labels):
        Z = mlab.bivariate_normal(
            X, Y, sigma[0], sigma[-1], mu[0], mu[1], sigma[1])
        ax_joint.contour(X, Y, Z, cmap=light_palette(
            color, as_cmap=True))
        Zx = mlab.normpdf(x, mu[0], sigma[0])
        ax_marg_x.plot(x, Zx, color=color, label=l)

        Zy = mlab.normpdf(y, mu[1], sigma[-1])
        ax_marg_y.plot(Zy, y, color=color)

        if l is not None:
            patches.append(mpatches.Patch(color=color, label=l))

    if len(patches) > 0:
        ax_joint.legend(handles=patches)

    if out_file is None:
        out_file = op.abspath('jointplot.pdf')

    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    return out_file
