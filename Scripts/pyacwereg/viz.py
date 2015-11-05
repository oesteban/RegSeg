#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2014-12-11 15:08:23
# @Last Modified by:   Oscar Esteban
# @Last Modified time: 2015-03-20 13:18:39
import os.path as op

def regseg_fig01(in_file, out_file):
    from matplotlib.path import Path
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt

    im = plt.imread(in_file)
    p1 = (880, 380)
    p2 = (265, 450)
    x = [p1[0], p2[0]]
    y = [p1[1], p2[1]]
    
    fig, ax = plt.subplots(figsize=(18,19))
    implot = ax.imshow(im)
    
    plt.xlim((180, 1250))
    plt.ylim((850, 0))
    
    
    ax.grid(True, color='b', alpha=0.2, linestyle='dashed', linewidth=5)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.xticks([600, 1150])
    plt.yticks([100, 650])
    
    a1 = ax.arrow(p1[0],p1[1],-70,40, head_width=20, head_length=50, fc='black', ec='black')
    a2 = ax.arrow(p2[0],p2[1],100,0, head_width=20, head_length=50, fc='black', ec='black')
    
    g1 = ax.arrow(p1[0],p1[1],10,-7, head_width=12, head_length=10, fc='r', ec='r')
    g2 = ax.arrow(p2[0],p2[1],50,0, head_width=12, head_length=10, fc='r', ec='r')
    
    ax.annotate(r"$\mathbf{v}_{1}$", xy=(p1[0]- 30, p1[1]+50), color='black', fontsize=60)
    ax.annotate(r"$\mathbf{v}_{2}$", xy=(p2[0]- 70, p2[1] + 10), color='black', fontsize=60)
    
    ax.annotate(r"$\hat{\mathbf{n}}_{1}$", xy=(p1[0]-100, p1[1]+ 110), color='black', fontsize=60)
    ax.annotate(r"$\hat{\mathbf{n}}_{2}$", xy=(p2[0]+100, p2[1]+80), color='black', fontsize=60)
    
    ax.annotate(r"$\bar{s}_{1}$", xy=(p1[0]+30, p1[1]), color='r', fontsize=60)
    ax.annotate(r"$\bar{s}_{2}$", xy=(p2[0]+10, p2[1]-25), color='r', fontsize=60)
    
    ax.annotate(r"$\Omega_{wm}$", xy=(370, 220), color='g', fontsize=90)
    ax.annotate(r"$\Omega_{gm}$", xy=(750, 220), color='darkblue', fontsize=90)
    ax.annotate(r"$\Omega_{bg}$", xy=(1100, 220), color='w', fontsize=90)
    
    ax.annotate(r"$\Gamma_0$", xy=(920, 610), color='g', fontsize=90)
    ax.annotate(r"$\Gamma_1$", xy=(1050, 570), color='darkblue', fontsize=90)
    
    bbox_props = dict(boxstyle="square", fc="w", ec="black", lw=2, alpha=0.7)
    ax.annotate(r"$\bar{s}_{1} = w_{0,1} \, \left[ \mathcal{D}^2_{wm} \left(M(\mathbf{v}_1)\right) - \mathcal{D}^2_{gm}\left(M(\mathbf{v}_1)\right)\right]\, \hat{\mathbf{n}}_1$", xy=(270, 800), color='black', fontsize=50, bbox=bbox_props)
    
    ax.scatter(x, y, color='black', edgecolor='w', s=150)
    
    ukx = [600, 600, 1150, 1150]
    uky = [100, 650, 100, 650]
    ax.scatter(ukx, uky, color='b', edgecolor='w', s=250, alpha=0.5)
    
    ax.annotate(r"$\mathbf{u}_{30}$", xy=(615, 120), color='b', fontsize=60, alpha=0.6)
    ax.annotate(r"$\mathbf{u}_{31}$", xy=(615, 670), color='b', fontsize=60, alpha=0.6)
    
    verts = [p1, (600, p1[1]), (600, 650)]
    codes = [ Path.MOVETO, Path.CURVE4, Path.CURVE4,]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', edgecolor='darkcyan', lw=2, alpha=0.4)
    ax.add_patch(patch)
    xs, ys = zip(*verts)
    ax.annotate(r"$\psi_{31}(\mathbf{v}_{1})}$", xy=(615, 590), color='darkcyan', fontsize=60, alpha=0.6)
    # ax.plot(xs, ys, 'x--', lw=2, color='b', ms=10)
    plt.savefig(out_file, format='pdf', bbox_inches='tight', pad_inches=0, dpi=300)

    return out_file


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
    g = sn.ConditionalJointGrid(f1name, f2name, df, hue='tissue',
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
                   huelabel='tissue', size=8, dpi=100, subsample=1.0,
                   xlims=None, ylims=None, out_file=None):
    import os.path as op
    import nibabel as nb
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sn
    import pandas as pd

    fa = nb.load(imageX).get_data().astype(np.float32).reshape(-1)
    md = nb.load(imageY).get_data().astype(np.float32).reshape(-1)
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

    if subsample < 1.0:
        df = df.sample(frac=subsample, replace=True)

    df = df[df.tissue != 'do-not-show']
    sn.set_context("talk", font_scale=1. + 0.1*size)

    g = sn.ConditionalJointGrid(xlabel, ylabel, df, hue=huelabel, xlim=xlims, ylim=ylims,
                     size=size, inline_labels=True)
    # g.plot_joint(plt.scatter)
    g.plot_joint(sn.kdeplot, linewidths=3.0)
    g.plot_marginals(sn.distplot)

    for ax in g.ax_joint:
        plt.setp(ax.get_yticklabels(), visible=False)
        plt.setp(ax.get_xticklabels(), visible=False)

    if out_file is not None:
        plt.savefig(out_file, dpi=dpi, bbox_inches='tight')
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


def slices_gridplot(in_files, view=['axial'], size=(5, 5), discard=3,
                    slices=None, label=None, out_file=None):
    from matplotlib import pyplot as plt
    import numpy as np

    view = np.atleast_1d(view).tolist()
    rows = len(view)

    fileslist = []
    if slices is None:
        cols = 0
        for v in view:
            filtlist = [f for f in in_files if v in f]
            if cols == 0:
                viewlist = filtlist[discard:-discard]
                cols = len(viewlist)
            else:
                viewlist = filtlist[discard:cols + discard]
            fileslist.append(viewlist)
    else:
        for v in view:
            cview = [f for f in in_files for s in slices
                     if '%s%04d' % (v, s) in f]
            fileslist.append(cview)

        rows, cols = np.shape(fileslist)

    fig, axes = plt.subplots(rows, cols,
                             figsize=(size[0] * cols, size[1] * rows),
                             subplot_kw={'xticks': [], 'yticks': []})
    fig.subplots_adjust(hspace=0.1, wspace=0.05)

    axes = np.atleast_2d(axes)
    for r in range(rows):
        for c in range(cols):
            ax = axes[r, c]
            im = plt.imread(fileslist[r][c])
            ax.imshow(im)

    if label is not None:
        label = np.atleast_1d(label).tolist()
        for i, l in enumerate(label):
            axes[i][0].set_ylabel(l, fontsize=40)

    if out_file is None:
        out_file = op.abspath('slices_gridplot.pdf')

    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    return out_file


def phantom_errors(in_csv, size=(80, 30), out_file=None):
    import pandas as pd
    import seaborn as sn
    import matplotlib.pyplot as plt
    sn.set_context("poster", font_scale=4.5)
    df = pd.read_csv(in_csv).drop_duplicates(
        subset=['surf_id', 'repetition', 'surfdist_avg', 'model_type'])
    del df['Unnamed: 0']

    df.surf_id[df.surf_id == 0] = 'internal'
    df.surf_id[df.surf_id == 1] = 'external'
    mtypes = df.model_type.unique()
    cols = len(mtypes) * 2
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=size)

    lodf = df[df.resolution == 'lo']
    plot0 = sn.violinplot(x='model_type', y='surfdist_avg', hue='surf_id',
                          hue_order=['internal', 'external'], inner='quartile',
                          data=lodf, scale_hue=.9, width=.9, ax=ax1)
    plot0.set_xlabel('Model Type')
    plot0.set_ylabel('Averaged error of surfaces (mm)')
    plot0.set_ylim([0.0, 2.25])
    plot0.set_title(
        r'Registration error @ $%.1f \times %.1f \times %.1fmm^3$' %
        tuple([2.0] * 3))

    l = plot0.axhline(y=2.0, lw=15, xmin=0.07, xmax=0.93,
                      color='gray', alpha=.4)
    plot0.annotate(
        "voxel size", xy=(2.5, 2.0), xytext=(100, -150),
        xycoords='data', textcoords='offset points', va='center',
        color='w', fontsize=80,
        bbox=dict(boxstyle='round', fc='gray', ec='none', color='w'),
        arrowprops=dict(arrowstyle='wedge,tail_width=.7',
                        fc='gray', ec='none',
                        ))

    frame = plot0.legend(loc=2, fancybox=True).get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('white')

    hidf = df[df.resolution == 'hi']
    plot1 = sn.violinplot(x='model_type', y='surfdist_avg', hue='surf_id',
                          hue_order=['internal', 'external'], inner='quartile',
                          data=hidf, scale_hue=.9, width=.9, ax=ax2)
    plot1.set_xlabel('Model Type')
    plot1.set_ylabel('')
    plot1.set_ylim([0.0, 1.15])
    plot1.set_title(
        r'Registration error @ $%.1f \times %.1f \times %.1fmm^3$' %
        tuple([1.0] * 3))
    l = plot1.axhline(y=1.0, lw=15, xmin=0.07, xmax=0.93,
                      color='gray', alpha=.4)
    plot1.annotate(
        "voxel size", xy=(0.1, 1.0), xytext=(-200, -120),
        xycoords='data', textcoords='offset points', va='center',
        color='w', fontsize=80,
        bbox=dict(boxstyle='round', fc='gray', ec='none', color='w'),
        arrowprops=dict(arrowstyle='wedge,tail_width=.7',
                        fc='gray', ec='none',
                        ))
    plot1.legend_.remove()

    sn.despine(left=True, bottom=True)

    if out_file is None:
        out_file = op.abspath('phantom_violinplots.pdf')

    plt.savefig(out_file, dpi=320, bbox_inches='tight')

#    lg = plot0.add_legend(title='Surface')
    return fig


def realdata_errors(in_csv, size=(80, 25), out_file=None,
                    columns=['surf_dist_1', 'surf_dist_3', 'surf_dist_5']):
    import pandas as pd
    import seaborn as sn
    import matplotlib.pyplot as plt
    import numpy as np

    def _extract_method(mdf, method_name, columns):
        rsd = mdf[columns].values.T
        values = rsd.reshape(-1).tolist()
        numvals = rsd.shape[1]
        totalvals = len(values)
        surflist = ([0] * numvals + [1] * numvals +
                    [2] * numvals + [3] * totalvals)

        data = {'method': [method_name] * (totalvals * 2),
                'surf': surflist,
                'error': values + values}
        return pd.DataFrame(data)

    fig = plt.figure(num=None, figsize=size)

    sn.set_context("poster", font_scale=5)
    df = pd.read_csv(in_csv).drop_duplicates(subset=['method', 'subject_id'])
    del df['Unnamed: 0']

    regsegdf = df[df.method == 'REGSEG'].reset_index(drop=True)
    t2bdf = df[df.method == 'T2B'].reset_index(drop=True)
    df1 = _extract_method(regsegdf, 'regseg', columns)
    df2 = _extract_method(t2bdf, 'T2B', columns)

    plot0 = sn.violinplot(x='surf', y='error', hue='method', size=size,
                          inner='quartile', linewidth=8,  # bw=0.2,
                          data=pd.concat([df1, df2]), scale_hue=.9, width=.7)

    ymax = plot0.get_ylim()[1]
    plot0.set_ylim([0.0, 4.0])

    plot0.set_xticklabels([r'$\Gamma_{VdGM}$', r'$\Gamma_{WM}$',
                           r'$\Gamma_{pial}$', 'Aggregated'])
    plot0.set_xlabel('Surface')
    plot0.set_ylabel('Surface warping index ($sWI$, mm)')
    l = plot0.axhline(y=1.25, lw=15, xmin=0.07, xmax=0.93,
                      color='gray', alpha=.4)
    plot0.annotate(
        "voxel size", xy=(-0.15, 1.25), xytext=(-250, 105),
        xycoords='data', textcoords='offset points', va='center',
        color='w', fontsize=80,
        bbox=dict(boxstyle='round', fc='gray', ec='none', color='w'),
        arrowprops=dict(arrowstyle='wedge,tail_width=.7',
                        fc='gray', ec='none',
                        ))

    plot0.set_title('Comparison in real data experiments')
    leg = plot0.legend(loc="best")
    sn.despine(left=True, bottom=True)

    if out_file is None:
        out_file = op.abspath('real_violinplots.pdf')
    plt.savefig(out_file, dpi=320, bbox_inches='tight')

#    lg = plot0.add_legend(title='Surface')
    return plot0


def phantom_boxplot(in_csv, resolution='lo',
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
                      size=15, kind='box', x_order=mtypes, legend=False,
                      n_boot=10e4, estimator=np.median)
    g.set_axis_labels('Type of model', 'Averaged error of surfaces (mm)')
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
    df = df.sort(columns=['subject_id', 'error']).drop_duplicates(
        subset=['subject_id', 'error', 'total'])
    subjects = df.subject_id.unique()

    series = []
    for s in subjects:
        ndf = df[df.subject_id == s]
        ndf.total /= ndf.total.min()
        ts = np.atleast_1d(ndf.total).tolist()
        series.append(ts)

    final = pd.DataFrame(np.array(series).T, index=ndf.error, columns=subjects)
    ax = final.plot(legend=False)
    ax.set_xlabel(r'Registration error ($\epsilon$)')
    ax.set_ylabel('Relative metric value')
    if out_file is not None:
        plt.savefig(out_file, dpi=300, bbox_inches='tight')

    return final
