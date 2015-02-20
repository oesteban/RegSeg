#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2014-12-11 15:08:23
# @Last Modified by:   Oscar Esteban
# @Last Modified time: 2015-02-20 15:29:00


def add_annotations(values, ax, level, nlevels, color, lastidx, units=''):
    labels = ('%.3g%s' % (values[0], units), '%.2g%s' % (values[1], units))
    xy_start = [0, values[0]]
    xy_end = [0, values[1]]

    xlength = -len(labels[0]) * 9.0
    if units != '':
        xlength -= 15
    xyt_start = (xlength, 0)
    xyt_end = xyt_start
    apos_start = (1.0, 0.5)
    apos_end = apos_start

    if (level == nlevels - 1):
        xy_start = [lastidx, values[0]]
        xy_end = [lastidx, values[1]]
        xyt_start = (15, 0)
        xyt_end = xyt_start
        apos_start = (0.0, 0.5)
        apos_end = apos_start

    if (level != nlevels - 1) and (level != 0):
        xy_start = [0, values[0]]
        xy_end = [lastidx, values[1]]
        xyt_start = (10, -25)
        xyt_end = (-len(labels[1]) * 9.5, 25)
        apos_start = (0.0, 1.0)
        apos_end = (1.0, 0.0)

    ax.annotate(
        labels[0], xy=tuple(xy_start), xytext=xyt_start,
        textcoords='offset points', va='center',
        color='w', fontsize='xx-small',
        bbox=dict(boxstyle='round', fc=color, ec='none', color='w'),
        arrowprops=dict(arrowstyle='wedge,tail_width=0.6',
                        fc=color, ec='none', relpos=apos_start,
                        ))

    ax.annotate(
        labels[1], xy=tuple(xy_end), xytext=xyt_end,
        textcoords='offset points', va='center',
        color='w', fontsize='xx-small',
        bbox=dict(boxstyle='round', fc=color, ec='none', color='w'),
        arrowprops=dict(arrowstyle='wedge,tail_width=0.6',
                        fc=color, ec='none', relpos=apos_end,
                        ))


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

    innersp = 0.08
    innergs0 = mgs.GridSpecFromSubplotSpec(
        3, len(levels), subplot_spec=outergs[0], hspace=0.6, wspace=innersp)

    ax3lims = (0, 0)
    for i, l in enumerate(levels):
        ax1 = plt.Subplot(fig, innergs0[0, i])
        ax2 = plt.Subplot(fig, innergs0[1, i])
        ax3 = plt.Subplot(fig, innergs0[2, i])

        fig.add_subplot(ax1)
        fig.add_subplot(ax2)
        fig.add_subplot(ax3)

        ldf = df[df.level == i]
        plt.sca(ax1)

        axplots = []
        axplots.append(ldf.E_t.plot(label=r'$E_T(t)$'))
        axplots.append(ldf.step.plot(
            secondary_y=True, style='k--', linewidth=0.7,
            label=r'$\delta(t)$', yticks=[]))

        its = ldf.iteration[ldf['desc_update'] != 0]
        if not its.empty:
            dupdt = ldf.desc_update[its] * ldf.E_t[its]
            ax1.plot(its, dupdt, marker='^', linestyle='',
                     fillstyle='full', color=palette[0])

        vals = (ldf.E_t[ldf.index[0]], ldf.E_t[ldf.index[-1]])
        add_annotations(vals, ax1, l, len(levels), palette[0], ldf.index[-1])

        if levels_df is not None:
            if not levels_df.empty:
                th_e = levels_df['total'][i]
                line, = ax1.plot([th_e] * (len(ldf['E_t'])))
                c = plt.getp(line, 'color')
                xy_start[1] = th_e
                ax1.annotate(
                    '%.3g' % th_e, xy=xy_end, xytext=xyt_end,
                    textcoords='offset points', va='center',
                    color='w', fontsize='xx-small',
                    bbox=dict(boxstyle='round', fc=c, ec='none', color='w'),
                    arrowprops=dict(arrowstyle='wedge,tail_width=0.6',
                                    fc=c, ec='none', relpos=apos_end,
                                    ))

        ax1.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax1.transAxes, fontdict=fd)
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.set_xlabel('')

        plt.sca(ax2)
        ldf.max_gk.plot(label=r'$\max_{k}\{\Vert \mathbf{g}_k(t) \Vert\}$')
        its = ldf.iteration[ldf['desc_update'] != 0]

        ldf.max_uk.plot(
            label=r'$\max_{k}\{\Vert \mathbf{u}_k(t) \Vert\}$', style='--')
        if not its.empty:
            dupdt = ldf.desc_update[its] * ldf.max_uk[its]
            ax2.plot(its, dupdt, marker='^', linestyle='',
                     fillstyle='full', color=palette[1])

        # vals = (ldf.max_gk[ldf.index[0]], ldf.max_gk[ldf.index[-1]])
        # add_annotations(vals, ax2, l, len(levels), palette[0], ldf.index[-1])

        vals = (ldf.max_uk[ldf.index[0]], ldf.max_uk[ldf.index[-1]])
        add_annotations(
            vals, ax2, l, len(levels), palette[1], ldf.index[-1], units='mm')

        ax2.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax2.transAxes, fontdict=fd)
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_xlabel('')

        if i == 0:
            hs = []
            ls = []
            for a in axplots:
                hlist, llist = a.get_legend_handles_labels()
                for h, l in zip(hlist, llist):
                    hs.append(h)
                    ls.append(l)
            bbox = (
                0.45, -.25, 0.75 * len(levels),
                1.0)
            lgd1 = ax1.legend(hs, ls, loc=3, ncol=len(ls), mode='expand',
                              bbox_to_anchor=bbox)
            lgd1.get_frame().set_facecolor('w')

            h, _ = ax2.get_legend_handles_labels()
            lgd2 = ax2.legend(loc=3, bbox_to_anchor=bbox, ncol=len(h),
                              mode='expand')
            lgd2.get_frame().set_facecolor('w')

        plt.sca(ax3)
        c = palette[0]
        c1 = palette[1]
        ldf.plot(x='iteration', y='g_50')
        ldf.plot(x='iteration', y='g_05', color=c1, style='--', lw=.7)
        ldf.plot(x='iteration', y='g_95', color=c1, style='--', lw=.7)
        ax3.fill_between(
            ldf.iteration, ldf['g_25'], ldf['g_75'], alpha=0.2, color=c)
        ax3.fill_between(
            ldf.iteration, ldf['g_05'], ldf['g_95'], alpha=0.05, color=c1)
        ax3.set_xlabel('iteration')

        if i == 0:
            ax3lims = ax3.get_ylim()
        else:
            ax3.set_ylim(ax3lims)

        if l == levels[-1]:
            ax3.yaxis.tick_right()
        if i == 0:
            ax3.set_ylabel(r'$g_i$ (mm)')

        if i > 0 and l != levels[-1]:
            ax3.set_yticklabels([])
            ax3.set_ylabel('')

        ax3.text(lposx, lposy, 'Level %d' %
                 (i + 1), ha=lha, transform=ax3.transAxes, fontdict=fd)

    innergs1 = mgs.GridSpecFromSubplotSpec(
        2, len(levels), subplot_spec=outergs[1], hspace=0.3, wspace=0.05)

    ax1lims = (0, 0)
    ax2lims = (0, 0)
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

        if l == levels[0]:
            ax1lims = ax1.get_ylim()
            ax2lims = ax2.get_ylim()
        else:
            ax1.set_ylim(ax1lims)
            ax2.set_ylim(ax2lims)

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
    fig.text(gps[2][0], gps[1][0] - 0.28,
             'Update evolution', size='medium')
    fig.text(gps[2][0], gps[0][0] + 0.2,
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


def jointplot_real(imageX, imageY, segmentation, mask,
                   labels=None, xlabel='X', ylabel='Y',
                   huelabel='tissue',
                   xlims=None, ylims=None, out_file=None):
    import os.path as op
    import nibabel as nb
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sn
    import pandas as pd

    fa = nb.load(imageX).get_data().reshape(-1)
    md = nb.load(imageY).get_data().reshape(-1)
    seg = nb.load(segmentation).get_data().astype(np.uint8).reshape(-1)
    msk = nb.load(mask).get_data().astype(np.uint8).reshape(-1)

    if xlims is None:
        xlims = (fa.min(), fa.max())

    if ylims is None:
        ylims = (md.min(), md.max())

    md[md < ylims[0]] = np.nan
    md[md > ylims[1]] = np.nan

    if labels is None:
        labels = ['Label%02d' % i for i in range(seg.max() + 1)]

    df_pieces = []
    for i, l in zip(np.unique(seg), labels):
        idxs = np.where((seg == i) & (msk > 0))
        fa_f = fa[idxs]
        md_f = md[idxs]
        d = {xlabel: fa_f,
             ylabel: md_f,
             huelabel: [l] * len(fa_f)}
        df_pieces.append(pd.DataFrame(d))

    df = pd.concat(df_pieces)
    df = df[df.tissue != 'do-not-show']
    sn.set_context("talk", font_scale=1.8)

    g = sn.JointGrid(xlabel, ylabel, df, hue=huelabel, xlim=xlims, ylim=ylims,
                     size=20, inline_labels=True)
    # g.plot_joint(plt.scatter)
    g.plot_joint(sn.kdeplot, linewidths=3.0)
    g.plot_marginals(sn.distplot)

    plt.setp(g.ax_joint.get_yticklabels(), visible=False)
    plt.setp(g.ax_joint.get_xticklabels(), visible=False)

    if out_file is not None:
        plt.savefig(out_file, dpi=300, bbox_inches='tight')
    return g


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


def slices_gridplot(in_files, view=['axial'], size=(3, 3), discard=3,
                    out_file=None):
    from matplotlib import pyplot as plt
    import numpy as np

    view = np.atleast_1d(view).tolist()
    rows = len(view)
    cols = 0

    fileslist = []
    for v in view:
        filtlist = [f for f in in_files if v in f]
        if cols == 0:
            viewlist = filtlist[discard:-discard]
            cols = len(viewlist)
        else:
            viewlist = filtlist[discard:cols + discard]
        fileslist.append(viewlist)

    fig, axes = plt.subplots(rows, cols,
                             figsize=(size[0] * cols, size[1] * rows),
                             subplot_kw={'xticks': [], 'yticks': []})
    fig.subplots_adjust(hspace=0.1, wspace=0.05)

    for r in range(rows):
        for c in range(cols):
            axes[r, c].imshow(plt.imread(fileslist[r][c]))

    if out_file is None:
        out_file = op.abspath('slices_gridplot.pdf')

    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    return out_file


def phantom_errors(in_csv, resolution='lo',
                   mtypes=['ball', 'box', 'L', 'gyrus']):
    import pandas as pd
    import seaborn as sn
    import matplotlib

    df = pd.read_csv(in_csv)
    df.surf_id[df.surf_id == 0] = 'internal'
    df.surf_id[df.surf_id == 1] = 'external'

    filtdf = df[df.resolution == resolution]

    sn.set(style="whitegrid")
    sn.set_context("poster", font_scale=1.5)
    g = sn.factorplot('model_type', 'surfdist_avg', 'surf_id', filtdf,
                      size=15, kind='box', x_order=mtypes, legend=False)
    g.set_axis_labels('Type of model', 'Averaged error of surfaces (mm)')
    lg = g.add_legend(title='Surface')

    mytitle = (r'Registration error @ $%.1f \times %.1f \times %.1fmm^3$' %
               tuple([1.0 if resolution == 'hi' else 2.0] * 3))
    g.fig.suptitle(mytitle, y=1.05, fontsize=30)

    return g


def metric_map_plot(in_file, out_file=None):
    import pandas as pd
    import numpy as np
    import seaborn as sn
    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    sn.set_context("talk", font_scale=1.25)
    plt.rcParams['font.size'] = 24

    df = pd.read_csv(in_file)
    df.error.fillna(0.0, inplace=True)
    df = df.sort(columns=['subject_id', 'error'])
    subjects = df.subject_id.unique()

    series = []
    for s in subjects:
        ndf = df[df.subject_id == s]
        ndf.total /= ndf.total.min()
        ts = ndf.total
        series.append(ts)

    final = pd.DataFrame(np.array(series).T, index=ndf.error, columns=subjects)
    ax = final.plot(legend=False)
    ax.set_xlabel(r'Registration error ($\epsilon$)')
    ax.set_ylabel('Relative metric value')
    if out_file is not None:
        plt.savefig(out_file, dpi=300, bbox_inches='tight')

    return final
