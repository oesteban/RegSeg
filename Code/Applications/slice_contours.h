// --------------------------------------------------------------------------------------
// File:          slice_contours.h
// Date:          Nov 18, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of ACWE-Reg
//
// ACWE-Reg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWE-Reg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWE-Reg.  If not, see <http://www.gnu.org/licenses/>.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef SLICE_CONTOURS_H_
#define SLICE_CONTOURS_H_


#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include <vector>

#include <vtkVersion.h>
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

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkImageToVTKImageFilter.h>


#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkWindowToImageFilter.h>

#if VTK_MAJOR_VERSION == 6
// https://gist.github.com/certik/5687727
// These 2 lines are needed due to:
// http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Factories_now_require_defines
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL);
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#else
#include <vtkGraphicsFactory.h>
#include <vtkImagingFactory.h>
#endif


namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

static const unsigned int DIMENSION = 3;

typedef itk::Image<float, DIMENSION>                         ImageType;
typedef typename ImageType::Pointer                          ImagePointer;
typedef itk::ImageFileReader<ImageType>                      ReaderType;
typedef typename ReaderType::Pointer                         ReaderPointer;
typedef itk::ImageToVTKImageFilter<ImageType>                ConnectorType;
typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleFilterType;

int main(int argc, char *argv[]);

#endif /* SLICE_CONTOURS_H_ */
