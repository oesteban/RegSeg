/* --------------------------------------------------------------------------------------
 * File:    regseg.h
 * Date:    02/10/2012
 * Author:  Oscar Esteban oesteban@die.upm.es
 * Version: 0.1
 * License: BSD
 * --------------------------------------------------------------------------------------

 Copyright (c) 2012, Oscar Esteban - oesteban@die.upm.es
 with Biomedical Image Technology, UPM (BIT-UPM) and
 Signal Processing Laboratory 5, EPFL (LTS5-EPFL).
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the names of the BIT-UPM and the LTS5-EPFL, nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#ifndef REGSEG_H_
#define REGSEG_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#include <itkVector.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkQuadEdgeMesh.h>
#include <itkVTKPolyDataReader.h>
#include <itkVTKPolyDataWriter.h>
#include <itkVectorImageToImageAdaptor.h>
#include <itkComposeImageFilter.h>
#include <itkVectorResampleImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkDisplacementFieldTransform.h>
#include <itkResampleImageFilter.h>
#include "MahalanobisFunctional.h"
#include "SpectralGradientDescentOptimizer.h"
#include "SpectralADMMOptimizer.h"
//#include "GradientDescentFunctionalOptimizer.h"
//#include "ALOptimizer.h"

using namespace rstk;


namespace bpo = boost::program_options;

typedef itk::Image<float, 3u>                                ChannelType;
typedef itk::Vector<float, 2u>                               VectorPixelType;
typedef itk::Image<VectorPixelType, 3u>                      ImageType;
typedef itk::ComposeImageFilter< ChannelType,ImageType >     InputToVectorFilterType;

typedef MahalanobisFunctional<ImageType>                      FunctionalType;
typedef FunctionalType::ContourType                           ContourType;
typedef ContourType::Pointer                                 ContourDisplacementFieldPointer;
typedef FunctionalType::MeanType                              MeanType;
typedef FunctionalType::CovarianceType                        CovarianceType;
typedef FunctionalType::FieldType                             DeformationFieldType;

typedef SpectralGradientDescentOptimizer< FunctionalType >   Optimizer;
//typedef SpectralADMMOptimizer< FunctionalType >              Optimizer;

typedef typename Optimizer::Pointer                          OptimizerPointer;

typedef itk::VTKPolyDataReader< ContourType >     ReaderType;
typedef itk::VTKPolyDataWriter< ContourType >     WriterType;
typedef itk::ImageFileReader<ChannelType>                      ImageReader;
typedef itk::ImageFileWriter<ChannelType>                      ImageWriter;
typedef itk::ImageFileWriter<DeformationFieldType>           DeformationWriter;

typedef itk::VectorImageToImageAdaptor<double,3u>            VectorToImage;
typedef itk::VectorResampleImageFilter
		<DeformationFieldType,DeformationFieldType,double>   DisplacementResamplerType;
typedef itk::BSplineInterpolateImageFunction
		                <DeformationFieldType>               InterpolatorFunction;
typedef itk::DisplacementFieldTransform<float, 3u>           TransformType;

typedef itk::ResampleImageFilter<ChannelType,ChannelType,float>    ResamplerType;



int main(int argc, char *argv[]);

#endif /* REGSEG_H_ */
