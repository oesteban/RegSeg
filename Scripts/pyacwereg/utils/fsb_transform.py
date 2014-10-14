import numpy as np
import nibabel as nib
from tvtk.api import tvtk
import datetime as dt
import os
import argparse
import sys
import subprocess as sp


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


def getfname(path):
    import os.path as op
    mname, mext = op.splitext(op.basename(path))
    if mext == ".gz":
        mname, _ = op.splitext(mname)
    return mname


def mri_info(fname, argument):
    import subprocess as sp
    import numpy as np

    cmd_info = "mri_info --%s %s" % (argument, fname)
    proc = sp.Popen(cmd_info, stdout=sp.PIPE, shell=True)
    data = bytearray(proc.stdout.read())
    mstring = np.fromstring(data.decode("utf-8"), sep='\n')
    result = np.reshape(mstring, (4, -1))
    return result


def transformVolume(mov_path, tar_path, tfmfile, invert=False, outname=''):
    from nipy.algorithms.resample import resample
    from nipy import load_image
    from nipy import save_image
    from nipy.core.api import Image
    import os
    import os.path as op

    # Load Images
    movimg = load_image(mov_path)
    tarimg = load_image(tar_path)

    # Load transform components
    mov_ras2vox_tkr = mri_info(mov_path, "ras2vox-tkr")
    tar_vox2ras_tkr = mri_info(tar_path, "vox2ras-tkr")
    reg = np.genfromtxt(tfmfile, skiprows=4, comments='r')

    # Compute transformation
    target_to_tkr = np.dot(
        tar_vox2ras_tkr, np.linalg.inv(nib.load(tar_path).get_affine()))
    tkr_to_moving = np.dot(nib.load(mov_path).get_affine(), mov_ras2vox_tkr)
    xfm = np.dot(tkr_to_moving, np.dot(np.linalg.inv(reg), target_to_tkr))

    # Resample T1 in dMRI space
    if not invert:
        result = resample(
            movimg, tarimg.coordmap, np.linalg.inv(xfm), np.shape(tarimg.get_data()))
        if outname == '':
            outname = op.join(op.dirname(
                mov_path), '%s-to-%s.nii.gz' % (getfname(mov_path), getfname(tar_path)))
    else:
        result = resample(
            tarimg, movimg.coordmap, xfm, np.shape(movimg.get_data()))
        if outname == '':
            outname = op.join(op.dirname(
                tar_path), '%s-to-%s.nii.gz' % (getfname(tar_path), getfname(mov_path)))

    d = save_image(result, outname)

    return outname


def MRIPreTess(in_file, in_norm, label=1, out_file=""):
    import subprocess as sp
    import os
    import os.path as op

    if out_file == "":
        fname, ext = op.splitext(op.basename(in_file))
        if ext == ".gz":
            fname, _ = op.splitext(fname)

        out_file = op.join(op.abspath(os.getcwd()), "%s_pretess.mgz" % fname)

    cmd = "mri_pretess %s %d %s %s" % (in_file, label, in_norm, out_file)
    proc = sp.check_call(cmd, shell=True)

    return out_file


def transformSurface(origin, mov_path, tar_path, tfmfile, out_file="", use_paraview=False, write_vtp=False, del_vtk=False):
    import numpy as np
    import nibabel as nib
    import os.path as op
    from tvtk.api import tvtk
    import subprocess as sp

    def mri_info(fname, argument):
        import subprocess as sp
        import numpy as np
        cmd_info = "mri_info --%s %s" % (argument, fname)
        proc = sp.Popen(cmd_info, stdout=sp.PIPE, shell=True)
        data = bytearray(proc.stdout.read())
        result = np.reshape(
            np.fromstring(data.decode("utf-8"), sep='\n'), (4, -1))
        return result

    # Load transform components
    mov_ras2vox_tkr = mri_info(mov_path, "ras2vox-tkr")
    tar_vox2ras_tkr = mri_info(tar_path, "vox2ras-tkr")
    ac = nib.load(tar_path).get_affine()
    if op.exists(tfmfile):
        reg = np.genfromtxt(tfmfile, skiprows=4, comments='r')
    else:
        reg = np.identity(4)

    # Compute transformation
    target_to_tkr = np.dot(tar_vox2ras_tkr, np.linalg.inv(ac))
    tkr_to_moving = np.dot(nib.load(mov_path).get_affine(), mov_ras2vox_tkr)
    xfm = np.dot(tkr_to_moving, np.dot(np.linalg.inv(reg), target_to_tkr))
    #xfm = np.dot( tkr_to_moving, np.linalg.inv(reg) )

    ofname, origin_ext = op.splitext(origin)

    if not origin_ext == 'vtk':
        # Call mris_convert
        new_origin = "%s_tmp.vtk" % origin
        sp.check_call('mris_convert %s %s' % (origin, new_origin), shell=True)
        origin = new_origin

        r = tvtk.PolyDataReader(file_name=origin)
        vtk = r.output
        r.update()

        nO = tkr_to_moving[0:3, 3]
        for i, point in enumerate(vtk.points):
            vtk.points[i] = point + nO

        fixed_center = "%s_fstkr.vtk" % ofname

        w = tvtk.PolyDataWriter(file_name=fixed_center, input=vtk)
        w.update()

        #savevtk( fsvtk, origin )
        # del_vtk==True

    else:
        r = tvtk.PolyDataReader(file_name=origin)
        vtk = r.output
        r.update()

    if use_paraview:
        M = xfm[0:3, 0:3]
        O = xfm[0:3, 3]

        phy_size = np.dot(
            ac[0:3, 0:3], np.array(nib.load(tar_path).get_data().shape))
        # print phy_size
        # M[0,0:3]*=-1
        # M[1,1]*=-1
        O[0] -= phy_size[0] * 0.5

    else:
        M = xfm[0:3, 0:3]
        O = xfm[0:3, 3]

    for i, point in enumerate(vtk.points):
        vtk.points[i] = np.dot(M, point) + O

    if write_vtp:
        writer = tvtk.XMLPolyDataWriter(
            file_name='%s.vtp' % output_prefix, input=vtk)
        writer.update()

    if out_file == '':
        ofname, origin_ext = op.splitext(origin)
        out_file = '%s_tfm.vtk' % ofname

    #savevtk( vtk, outfname )
    w = tvtk.PolyDataWriter(file_name=out_file, input=vtk)
    w.update()

    if del_vtk:
        os.remove(origin)

    return out_file


def main():
    import sys
    parser = argparse.ArgumentParser(
        description="Transform freesurfer's surface to a targ space using an affine transformation")
    parser.add_argument('--targ', dest='target', type=lambda x: is_valid_file(
        parser, x), help='target volume (e.g. b=0 volume in dMRI)', required=True)
    parser.add_argument('--orig', dest='origin', type=lambda x: is_valid_file(parser, x),
                        help='original surface (e.g. surf/lh.pial outcome from freesurfer)', required=True)
    parser.add_argument('--o', dest='output_prefix', type=str,
                        help='output prefix for filename (will be an XML Vtk file)')
    parser.add_argument('--tfm', dest='transform', type=lambda x: is_valid_file(
        parser, x), help='transform file (e.g. register.dat from bbregister)', required=True)
    parser.add_argument('--paraview', dest='use_paraview', action='store_true',
                        default=False, help='True if orientation should be fixed for Paraview rendering')
    args = parser.parse_args()

    if args.output_prefix == None:
        basename, ext = os.path.splitext(args.target)
        args.output_prefix = [basename]
        if ext == '.gz':
            args.output_prefix, _ = os.path.splitext(basename)
        args.output_prefix += '_%s' % os.path.basename(args.origin)

    transform(args.origin, args.target, args.transform,
              args.output_prefix, args.use_paraview)

    return 0


def savevtk(vtk, fname, origin=(0.0, 0.0, 0.0)):
    ncells = vtk.polys.number_of_cells
    npoints = len(vtk.points)
    nvalues = len(vtk.polys.data)
    cells = np.reshape(vtk.polys.data, (ncells, int(nvalues / ncells)))
    now = dt.datetime.now()
    points = vtk.points + np.array(origin)

    if os.path.exists(fname):
        os.remove(fname)

    try:
        with open(fname, 'a+') as f:
            f.write('# vtk DataFile Version 2.0\n')  # Write identification
            f.write('%s, generated on %s\n' % (os.path.basename(fname), now.strftime(
                "%Y-%m-%d %H:%M")))  # Write comment (up to 256 chars)
            f.write('ASCII\n')
            f.write('DATASET POLYDATA\n')
            f.write('POINTS %d float\n' % npoints)
            for row in points:
                f.write('%.6f %.6f %.6f\n' % tuple(row))

            f.write('POLYGONS %d %d\n' % (ncells, nvalues))
            for row in cells:
                f.write('3 %d %d %d\n' % tuple(row[1:]))
    except Exception as err:
        print err
        return 0

    return 1


if __name__ == '__main__':
    sys.exit(main())


def disk_structure(size=(2, 2, 2)):
    totalSize = np.array(size) * 2 + 1
    struct = np.zeros(totalSize)
    coord = np.indices(totalSize).astype(np.float32)
    positions = coord.T - ((totalSize.astype(np.float32) - 1) / 2)
    terms = np.sum(
        ((positions ** 2) / (np.array(size) ** 2).astype(np.float32)).T, axis=0)
    struct[terms <= 1] = 1
    return struct.astype(np.bool)


def CombineSegs(in_file, labels):
    import nibabel as nib
    import numpy as np
    import os.path as op
    from scipy import ndimage

    img = nib.load(in_file)
    imgdata = img.get_data()

    outdata = np.zeros(shape=imgdata.shape)
    for l in labels:
        outdata[imgdata == l] = 1

    #disk = disk_structure( (2,2,2) )
    #outdata = ndimage.binary_opening( outdata, structure=disk ).astype(np.int)
    #outdata = ndimage.binary_closing( outdata, structure=disk ).astype(np.int)

    fname, ext = op.splitext(in_file)

    if ext == '.gz':
        fname, _ = op.splitext(fname)

    out_file = '%s_bin.nii.gz' % fname
    nib.save(
        nib.Nifti1Image(outdata, img.get_affine(), img.get_header()), out_file)
    return out_file

# 3

import os
import numpy as np
from tvtk.api import tvtk
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from pandas import *
import nibabel as nib
import matplotlib.cm as cm
import datetime as dt
import os

import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as niu
import nipype.interfaces.io as nio


def surfaceWorkflow(model_path, name="Surface"):
    pipeline = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_roi', 'surf_id', 'model_id']), name='inputnode')
    tess = pe.Node(fs.MRITessellate(
        label_value=1, use_real_RAS_coordinates=True), name='tessellate')
    smooth = pe.Node(
        fs.SmoothTessellation(disable_estimates=True), name='smooth')
    convert = pe.Node(fs.MRIsConvert(out_datatype='vtk'), name='tovtk')
    ds = pe.Node(nio.DataSink(base_directory=model_path), name='sinker')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['surface']), name='outputnode')

    pipeline.connect([
        (inputnode, tess, [('in_roi', 'in_file')]), (tess, smooth, [('surface', 'in_file')]), (smooth, convert, [('surface', 'in_file')]), (inputnode, ds, [
            ('model_id', 'container'), ('surf_id', 'surf')]), (convert, ds, [('converted', 'surf.@file')]), (ds, outputnode, [('out_file', 'surface')])
    ])
    return pipeline


def fixVtk(input_file, output_file, ref_nii):
    ref = nib.load(ref_nii)
    orig = np.array(ref.get_shape()) * 0.5
    ac = ref.get_affine()[0:3, 0:3]

    rotM = np.zeros(shape=ac.shape)
    rotM[0, 0] = -1.0
    rotM[1, 2] = 1.0
    rotM[2, 1] = -1.0

    R = np.dot(rotM, ac)

    with open(input_file, 'r') as f:
        with open(output_file, 'w') as w:
            npoints = 0
            pointid = -5

            for i, l in enumerate(f):
                if (i == 4):
                    s = l.split()
                    npoints = int(s[1])
                    fmt = np.dtype(s[2])
                elif(i > 4 and pointid < npoints):
                    vert = np.array([float(x) for x in l.split()])
                    vert = np.dot(vert, R) + orig
                    l = '%.9f  %.9f  %.9f\n' % tuple(vert)

                w.write(l)
                pointid = pointid + 1
    return True


def GenerateSurface(data, fname, use_smooth=True, spacing=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0)):
    grid = tvtk.ImageData(
        spacing=spacing, origin=origin, dimensions=data.shape)
    grid.point_data.scalars = data.ravel().astype(np.uint8)
    grid.point_data.scalars.name = 'scalars'
    iso = tvtk.ImageMarchingCubes(input=grid)
    iso.update()
    contour = iso.output
    if use_smooth:
        vtkmeshclean = tvtk.CleanPolyData(input=iso.output)
        vtkmeshclean.update()
        smooth = tvtk.SmoothPolyDataFilter(
            input=vtkmeshclean.output, number_of_iterations=20)
        smooth.update()
        vtkmeshclean2 = tvtk.CleanPolyData(input=smooth.output)
        vtkmeshclean2.update()
        contour = vtkmeshclean2.output
    contour.point_data.scalars = np.ones(np.shape(contour.points)[0])
    contour.point_data.scalars.name = 'scalars'
    # origin = tuple( (np.array(data.shape).astype(float)) * -0.5 ) # fix new
    # origin to have a nice view in tkmedit
    vtk2fsa(contour, '%s.asc' % os.path.splitext(fname)[0], origin=origin)
    #w = tvtk.PolyDataWriter(input=contour, file_name=fname )
    # w.write()
    savevtk(contour, fname)
    return contour


def vtk2fsa(vtk, fname, origin=(0.0, 0.0, 0.0)):
    ncells = vtk.polys.number_of_cells
    npoints = len(vtk.points)
    nvalues = len(vtk.polys.data)
    cells = np.reshape(vtk.polys.data, (ncells, int(nvalues / ncells)))
    now = dt.datetime.now()
    points = vtk.points + np.array(origin)

    if os.path.exists(fname):
        os.remove(fname)

    try:
        with open(fname, 'a+') as f:
            f.write('#!ascii freesurfer mesh generated from vtk - %s, generated on %s\n' %
                    (os.path.basename(fname), now.strftime("%Y-%m-%d %H:%M")))  # Write header
            # Write number of vertices and number of faces
            f.write('%d %d\n' % (npoints, ncells))

            for row in points:
                f.write('%.6f %.6f %.6f 0\n' % tuple(row))

            for row in cells:
                f.write('%d %d %d 0\n' % tuple(row[1:]))
    except Exception as err:
        print err
        return 0

    return 1


def FilterMesh(points, sl, zdim, fname=None):
    lower = sl <= zdim * 0.5
#    filtered = np.around( points[(points[:,2]>(sl-0.25)) & (points[:,2]<(sl+0.25))  ] )
    if lower:
        filtered = np.around(points[(points[:, 2] < (sl + 0.25))])
    else:
        filtered = np.around(points[(points[:, 2] > (sl - 0.25))])

#    filtered = points[(points[:,2]>(sl-0.25))]
#    filtered = np.array(list(set(tuple(p) for p in filtered)))
    filtered.view('f,f,f').sort(order=['f2'], axis=0)
    uniq = np.array(DataFrame(filtered).drop_duplicates(
        cols=[0, 1], take_last=lower).values)
    if not fname == None:
        mesh = tvtk.PolyData(points=uniq)
        mesh.point_data.scalars = np.zeros(np.shape(uniq)[0], 'f')
        mesh.point_data.scalars.name = 'scalars'
        w = tvtk.PolyDataWriter(input=mesh, file_name=fname)
        w.write()
    return uniq, lower


def ComputeContourFunction(image, surf, sl, axis=2, fig=None, color='r'):
    points = np.array(surf.points, dtype=np.float32)
    if axis == 1:
        points[:, [0, 1, 2]] = points[:, [0, 2, 1]]
    elif axis == 0:
        points[:, [0, 1, 2]] = points[:, [2, 1, 0]]

    dim = np.shape(np.swapaxes(image.get_data(), axis, 2))
    dim2d = (dim[0], dim[1])
    xi, yi = np.mgrid[0:dim2d[0], 0:dim2d[1]]
    ps, lower = FilterMesh(points, sl, dim[2], 'Model1/test_pandas.vtk')
    if lower:
        z_level = sl - 0.05
    else:
        z_level = sl + 0.05
    interpolated = griddata(
        ps[:, [0, 1]], ps[:, 2], (xi, yi), method='cubic', fill_value=z_level)
    return interpolated, z_level


def PlotContour(ref_path, surfs, ax=2, sl=50, colors=['b', 'r', 'g', 'y']):
    image = nib.load(ref_path)
    sliced = np.swapaxes(image.get_data(), ax, 2)[:, :, sl]
    cs = [plt.imshow(sliced, cmap=cm.Greys_r)]

    for i, s in enumerate(surfs):
        pts, z_level = ComputeContourFunction(image, s, sl, axis=ax)
        cs.append(
            plt.contour(pts, linewidths=2.5, colors=colors[i], levels=[z_level]))
    return cs
