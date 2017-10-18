// This file is part of RegSeg
//
// Copyright 2014-2017, Oscar Esteban <code@oscaresteban.es>
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.

#ifndef SLICE_CONTOURS_H_
#define SLICE_CONTOURS_H_

// See: http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Factories_now_require_defines
// See: https://gist.github.com/certik/5687727
#include <vtkVersion.h>

#if VTK_MAJOR_VERSION <= 5
	#include <vtkGraphicsFactory.h>
	#include <vtkImagingFactory.h>
#else
	#include <vtkAutoInit.h>
#endif


#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include <vector>

#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkImageActor.h>
#include <vtkImageProperty.h>
#include <vtkImageMapper3D.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageData.h>
#include <vtkSphereSource.h>
#include <vtkMetaImageWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkStripper.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkImageMapper.h>
#include <vtkImageSliceMapper.h>
#include <vtkImageSlice.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPNGWriter.h>
#include <vtkCamera.h>
#include <vtkCornerAnnotation.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkMatrix4x4.h>
#include <vtkLookupTable.h>
#include <vtkImageMapToColors.h>
#include <vtkImageResliceMapper.h>


#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkRescaleIntensityImageFilter.h>

#include <vtkImageMagnify.h>
#include <vtkImageViewer2.h>
#include <vtkImageReslice.h>

#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkWindowToImageFilter.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkAxesActor.h>
#include <vtkCubeAxesActor2D.h>

#include "itkImageToVTKImageFilter.h"

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

static const unsigned int DIMENSION = 3;

typedef itk::Image<float, DIMENSION>                         ImageType;
typedef typename ImageType::Pointer                          ImagePointer;
typedef itk::ImageFileReader<ImageType>                      ReaderType;
typedef typename ReaderType::Pointer                         ReaderPointer;
typedef itk::ImageToVTKImageFilter<ImageType>                ConnectorType;
typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleFilterType;

// Matrices for axial, coronal, sagittal, oblique view orientations
static const double axialElements[16] = {
         1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1 };

static const double coronalElements[16] = {
         1, 0, 0, 0,
         0, 0, 1, 0,
         0,-1, 0, 0,
         0, 0, 0, 1 };

static const double sagittalElements[16] = {
         0, 0,-1, 0,
         1, 0, 0, 0,
         0,-1, 0, 0,
         0, 0, 0, 1 };

int main(int argc, char *argv[]);

#endif /* SLICE_CONTOURS_H_ */
