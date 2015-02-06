#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2015-02-06 13:26:12
# @Last Modified by:   oesteban
# @Last Modified time: 2015-02-06 13:26:25


def laplacian_filter(in_file, in_mask=None, out_file=None):
    import numpy as np
    import nibabel as nb
    import os.path as op
    from math import pi
    from numpy.fft import fftn, ifftn

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_smooth.nii.gz' % fname)

    im = nb.load(in_file)
    data = im.get_data()

    if in_mask is not None:
        mask = nb.load(in_mask).get_data()
        mask[mask > 0] = 1.0
        mask[mask <= 0] = 0.0
        data *= mask

    dataft = fftn(data)
    x = np.linspace(0, 2 * pi, dataft.shape[0])[:, None, None]
    y = np.linspace(0, 2 * pi, dataft.shape[1])[None, :, None]
    z = np.linspace(0, 2 * pi, dataft.shape[2])[None, None, :]
    lapfilt = 2.0 * np.squeeze((np.cos(x) + np.cos(y) + np.cos(z))) - 5.0
    dataft *= lapfilt
    imfilt = ifftn(dataft)

    nb.Nifti1Image(imfilt.astype(np.float32), im.get_affine(),
                   im.get_header()).to_filename(out_file)
    return out_file


def rbf_approx(in_file, in_mask=None, out_file=None):
    import numpy as np
    import nibabel as nb
    import os.path as op
    from scipy.interpolate import Rbf

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_interp.nii.gz' % fname)

    im = nb.load(in_file)
    aff = im.get_affine()
    imdata = im.get_data()
    result = np.zeros_like(imdata)

    if in_mask is not None:
        mask = nb.load(in_mask).get_data()
        mask[mask > 0] = 1.0
        mask[mask <= 0] = 0.0
        imdata *= mask
    else:
        mask = np.zeros_like(imdata)
        mask[imdata != 0] = 1

#    for z in np.arange(imdata.shape[2]):
#        zslice = imdata[:, :, z]
#        ij = np.where(np.ones_like(zslice))
#        data = zslice[ij]
#
#        aff = aff[0:2, 0:2]
#        xy = aff.dot(np.array([col for col in ij], dtype=np.float32))
#        rbfi = Rbf(xy[0, :], xy[1, :], data, function='gaussian')
#
#        grid = np.where(np.ones_like(zslice))
#        xy2 = aff.dot(np.array([col for col in grid], dtype=np.float32))
#
#        di = rbfi(xy2[0, :], xy2[1, :])
#        result[:, :, z] = di.reshape(zslice.shape)

    ijk = np.where(np.ones_like(imdata))
    rndidx = np.random.choice(len(ijk[0]), size=400, replace=False)
    ijk = tuple([cs[rndidx] for cs in ijk])

    data = imdata[ijk]

    xyz = aff.dot(np.array([col for col in ijk] + [[1] * len(ijk[0])],
                           dtype=np.float32))
    rbfi = Rbf(xyz[0, :], xyz[1, :], xyz[2, :], data,
               smooth=0.05, epsilon=60.0, function='inverse')

    grid = np.where(np.ones_like(imdata))
    xyz2 = aff.dot(np.array([col for col in grid] + [[1] * len(grid[0])],
                            dtype=np.float32))
    di = rbfi(xyz2[0, :], xyz2[1, :], xyz2[2, :])
    result = di.reshape(imdata.shape)

    result *= mask
    nb.Nifti1Image(result, im.get_affine(),
                   im.get_header()).to_filename(out_file)

    return out_file


def deconv(in_file, in_mask=None, out_file=None):
    import numpy as np
    import nibabel as nb
    import os.path as op
    from math import pi, sqrt, exp
    from scipy.stats import multivariate_normal as mnormal
    from scipy.ndimage.filters import gaussian_filter
    import numpy.linalg as nla

    def deconvolve(data, psf):
        from numpy.fft import fftn, ifftn, fftshift, ifftshift
        dataft = fftshift(fftn(data))
        psfft = fftshift(fftn(psf))
        return fftshift(ifftn(ifftshift(dataft / psfft)))

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_deconv.nii.gz' % fname)

    im = nb.load(in_file)
    data = im.get_data()

    mask = np.ones_like(data)
    if in_mask is not None:
        mask = nb.load(in_mask).get_data()
        mask[mask > 0] = 1.0
        mask[mask <= 0] = 0.0

    ijk = np.where(np.ones_like(data))
    xyz = [tuple(v) for v in np.asarray(ijk, dtype=float).T]

    di = mnormal.pdf(xyz, mean=np.array(data.shape) * 0.0, cov=np.eye(3) * 20)
    imfilt = deconvolve(data, di.reshape(data.shape)) * mask
    data = gaussian_filter(imfilt, sigma=10.0)

    nb.Nifti1Image(data.astype(np.float32), im.get_affine(),
                   im.get_header()).to_filename(out_file)
    return out_file


def wavelets_denoise(in_file, in_mask=None, out_file=None):
    import numpy as np
    import nibabel as nb
    import os.path as op
    from math import pi
    from scipy.interpolate import Rbf
    from scipy.ndimage import median_filter
    import pywt as wt

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_wavelets.nii.gz' % fname)

    im = nb.load(in_file)
    aff = im.get_affine()
    imdata = im.get_data()

    datamax = imdata.max()
    datamin = imdata.min()

    imdata = 255. * (imdata - datamin) / (datamax - datamin)

    if in_mask is not None:
        mask = nb.load(in_mask).get_data()
        mask[mask > 0] = 1.0
        mask[mask <= 0] = 0.0
        imdata *= mask
    else:
        mask = np.zeros_like(imdata)
        mask[imdata != 0] = 1

    result = np.zeros_like(imdata)
    wavelet = wt.Wavelet('db10')
    thres = 100

    offset = (0 if imdata.shape[0] % 2 == 0 else 1,
              0 if imdata.shape[1] % 2 == 0 else 1)

    for z in np.arange(imdata.shape[2]):
        zslice = imdata[offset[0]:, offset[1]:, z]
        wcoeff = wt.wavedec2(zslice, wavelet)
        nwcoeff = map(lambda x: wt.thresholding.soft(x, thres), wcoeff)
        result[offset[0]:, offset[1]:, z] = wt.waverec2(nwcoeff, wavelet)

    m = np.median(result[offset[0]:, offset[1]:, :])
    result[offset[0]:, offset[1]:, :] = 2.0 * \
        (result[offset[0]:, offset[1]:, :] - m) / 255.

    result = median_filter(result, 10)
    result *= (pi / np.percentile(result, 99.99))
    # result = (result / 255.) * (datamax - datamin) - datamin
    nb.Nifti1Image(result, im.get_affine(),
                   im.get_header()).to_filename(out_file)

    return out_file
