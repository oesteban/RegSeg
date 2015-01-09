#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2014-12-11 15:08:23
# @Last Modified by:   oesteban
# @Last Modified time: 2015-01-09 02:35:28


def plot_report(df, levels_df=None, out_file=None):
    import nibabel as nb
    import numpy as np
    import pandas as pd
    import seaborn as sn
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as mgs
    from matplotlib import patches as mpatches
    from itertools import cycle

    plt.style.use('ggplot')
    sn.set_context("talk", font_scale=1.25)
    plt.rcParams['font.size'] = 24
    palette = plt.rcParams['axes.color_cycle']
    levels = pd.Series(df.level).unique()

    fd = {'size': 'medium', 'color': 'w', 'weight': 'bold'}
    lposx = .05
    lposy = .9
    lha = 'left'

    fig = plt.figure(figsize=(20, 15))
    outergs = mgs.GridSpec(1, 2, wspace=0.2, hspace=0.2)

    innergs0 = mgs.GridSpecFromSubplotSpec(
        3, len(levels), subplot_spec=outergs[0], hspace=0.3, wspace=0.08)

    for i, l in enumerate(levels):
        isRight = ((i % 2) != 0)

        ax1 = plt.Subplot(fig, innergs0[0, i])
        ax2 = plt.Subplot(fig, innergs0[1, i])
        ax3 = plt.Subplot(fig, innergs0[2, i])

        fig.add_subplot(ax1)
        fig.add_subplot(ax2)
        fig.add_subplot(ax3)

        ldf = df[df.level == i]
        plt.sca(ax1)
        ldf.E_t.plot(label=r'$E_T(t)$')

        th_e = np.squeeze(ldf.tail(1)['E_t'].values)
        xy = [0, th_e]
        xytext = (-80, 0)
        apos = (1.0, 0.5)

        if isRight:
            xy = [ldf.index[-1], th_e]
            xytext = (15, 0)
            apos = (0.0, 0.5)

        ax1.annotate(
            '%.3g' % th_e, xy=tuple(xy), xytext=xytext,
            textcoords='offset points', va='center',
            color='w', fontsize='xx-small',
            bbox=dict(boxstyle='round', fc=palette[0], ec='none', color='w'),
            arrowprops=dict(arrowstyle='wedge,tail_width=0.6',
                            fc=palette[0], ec='none', relpos=apos,
                            ))

        th_e = np.squeeze(ldf.head(1)['E_t'].values)
        xy[1] = th_e
        ax1.annotate(
            '%.3g' % th_e, xy=tuple(xy), xytext=xytext,
            textcoords='offset points', va='center',
            color='w', fontsize='xx-small',
            bbox=dict(boxstyle='round', fc=palette[0], ec='none', color='w'),
            arrowprops=dict(arrowstyle='wedge,tail_width=0.6',
                            fc=palette[0], ec='none', relpos=apos,
                            ))

        its = ldf.iteration[ldf['desc_update'] != 0]
        if not its.empty:
            dupdt = ldf.desc_update[its] * ldf.E_t[its]
            ax1.plot(its, dupdt, marker='^', linestyle='',
                     fillstyle='full', color=palette[0])

        if levels_df is not None:
            if not levels_df.empty:
                th_e = levels_df['total'][i]
                line, = ax1.plot([th_e] * (len(ldf['E_t'])))
                c = plt.getp(line, 'color')

                xy[1] = th_e

                ax1.annotate(
                    '%.3g' % th_e, xy=xy, xytext=xytext,
                    textcoords='offset points', va='center',
                    color='w', fontsize='xx-small',
                    bbox=dict(boxstyle='round', fc=c, ec='none', color='w'),
                    arrowprops=dict(arrowstyle='wedge,tail_width=0.6',
                                    fc=c, ec='none', relpos=apos,
                                    ))

        ldf.step.plot(secondary_y=True, style='k--', linewidth=0.7,
                      label=r'$\delta(t)$', yticks=[])

        ax1.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax1.transAxes, fontdict=fd)
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.set_xlabel('')

        plt.sca(ax2)
        ldf.max_gk.plot(label=r'$\max_{k}\{\Vert \mathbf{g}_k(t) \Vert\}$')

        its = ldf.iteration[ldf['desc_update'] != 0]
        if not its.empty:
            dupdt = ldf.desc_update[its] * ldf.max_gk[its]
            ax2.plot(its, dupdt, marker='^', linestyle='',
                     fillstyle='full', color=palette[0])

        ldf.max_uk.plot(secondary_y=True, yticks=[],
                        label=r'$\max_{k}\{\Vert \mathbf{u}_k(t) \Vert\}$')

        ax2.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax2.transAxes, fontdict=fd)
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_xlabel('')

        plt.sca(ax3)
        ldf.plot(x='iteration', y='g_50')
        ldf.plot(x='iteration', y='g_05', color=c, style='--', lw=.7)
        ldf.plot(x='iteration', y='g_95', color=c, style='--', lw=.7)
        ax3.fill_between(
            ldf.iteration, ldf['g_25'], ldf['g_75'], alpha=0.2, color=c)
        ax3.fill_between(
            ldf.iteration, ldf['g_05'], ldf['g_95'], alpha=0.05, color=c)
        ax3.set_xlabel('iteration')

        if (i % 2) != 0:
            ax3.yaxis.tick_right()
        else:
            ax3.set_ylabel(r'$g_i$ (mm)')

        ax3.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax3.transAxes, fontdict=fd)

    innergs1 = mgs.GridSpecFromSubplotSpec(
        2, len(levels), subplot_spec=outergs[1], hspace=0.3, wspace=0.05)
    for i, l in enumerate(levels):
        ldf = df[df.level == i]
        ax1 = plt.Subplot(fig, innergs1[0, i])
        ax2 = plt.Subplot(fig, innergs1[1, i])
        fig.add_subplot(ax1)
        fig.add_subplot(ax2)

        eid = 0
        erois = []

        while 'E_%02d' % eid in ldf:
            erois.append('E_%02d' % eid)
            eid += 1
        erois = erois[:-1]

        patches1 = []
        for n in range(len(erois)):
            patches1.append(
                mpatches.Patch(color=palette[n], label='$E_{%02d}$' % n))

        ldf.plot(ax=ax1, x='iteration', y=erois,
                 kind='area', stacked=True, yticks=[])
        ax1.text(.9, .05, 'Level %d' %
                 (i + 1), ha='right', transform=ax1.transAxes, fontdict=fd)

        plt.sca(ax2)
        ldf.plot(ax=ax2, x='iteration', y=['E_d', 'E_r'],
                 kind='area', stacked=True, yticks=[])
        plt.plot(ldf['E_t'])
        ax2.text(.9, .05, 'Level %d' %
                 (i + 1), ha='right', transform=ax2.transAxes, fontdict=fd)

    legend = ax1.legend(loc='center right', handles=patches1)
    frame = legend.get_frame()
    frame.set_facecolor('w')

    patches2 = []
    for n, txt in enumerate(['$E_{data}$', '$E_{reg}$']):
        patches2.append(mpatches.Patch(color=palette[n], label=txt))

    legend = ax2.legend(loc='center right', handles=patches2)
    frame = legend.get_frame()
    frame.set_facecolor('w')

    fig.suptitle('Convergence report', fontsize='x-large')

    gps = outergs.get_grid_positions(fig)
    fig.text(gps[2][0], gps[1][0] + 0.02,
             'Total energy evolution', size='medium')
    fig.text(gps[2][0], gps[1][0] - 0.26,
             'Update evolution', size='medium')
    fig.text(gps[2][0], gps[0][0] + 0.24,
             r'Gradient at vertices ($g_i$) distribution', size='medium')

    fig.text(gps[3][0] + 0.072, gps[1][0] + 0.02,
             'Region-wise evolution of energy', size='medium')
    fig.text(gps[3][0] + 0.072, gps[1][0] - 0.41,
             'Energy decomposition', size='medium')

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
                     xlim=f1lims, ylim=f2lims, size=20,
                     inline_labels=True)
    g.plot_joint(sn.kdeplot, linewidths=3.0)
    g.plot_marginals(sn.distplot)

    nlabels = len(g.ax_joint.yaxis.get_ticklabels())

    plt.setp(g.ax_joint.get_yticklabels(), visible=False)
    plt.setp(g.ax_joint.get_xticklabels(), visible=False)

    if out_file is None:
        out_file = op.abspath('jointplot.pdf')

    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    return g, out_file


def jointplot_gmm(locs, covs, labels=None, out_file=None,
                  xlims=None, ylims=None,
                  xname='Image A', yname='Image B',
                  size=20, ratio=5, space=.2):
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
    from scipy.stats import multivariate_normal

    if len(locs) != len(covs):
        raise RuntimeError('Mixture model does not contain elements')

    ncomp = np.shape(locs)[1]

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
    palette = color_palette("husl", n_colors=len(locs))

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

    X, Y = np.mgrid[
        xlims[0]:xlims[1] + xstep:xstep, ylims[0]:ylims[1] + ystep:ystep]

    pos = np.empty(X.shape + (2,))
    pos[:, :, 0] = X
    pos[:, :, 1] = Y
    y = Y[0]
    x = pos.T[0][0]

    if labels is None:
        labels = [None] * len(locs)

    patches = []
    for mu, cov, color, l in zip(locs, covs, palette, labels):
        cov = np.array(cov).reshape(ncomp, ncomp)
        mv = multivariate_normal(mu, cov)
        ZC = ax_joint.contour(X, Y, mv.pdf(pos),
                              cmap=light_palette(color, as_cmap=True))

        if l is not None:
            ax_joint.annotate(
                l, xy=(mu[0], mu[1]), xytext=(30, 20),
                textcoords='offset points', size=30, va='center',
                color='w',
                bbox=dict(boxstyle="round", fc=color, ec='none', alpha=0.9,
                          color='w'),
                arrowprops=dict(arrowstyle="wedge,tail_width=0.7",
                                fc=color, ec="none", alpha=0.6,
                                relpos=(0.2, 0.5),
                                )
            )

        Zx = mlab.normpdf(x, mu[0], np.sqrt(cov[0][0]))
        ax_marg_x.plot(x, Zx, color=color, label=l)

        Zy = mlab.normpdf(y, mu[1], np.sqrt(cov[1][1]))
        ax_marg_y.plot(Zy, y, color=color)

    if out_file is None:
        out_file = op.abspath('jointplot.pdf')

    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    return out_file
