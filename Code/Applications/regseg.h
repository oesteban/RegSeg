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

#ifndef REGSEG_H_
#define REGSEG_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>


#include <itkVector.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkThresholdImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkQuadEdgeMesh.h>
#include <itkVTKPolyDataReader.h>
#include "rstkVTKPolyDataWriter.h"
#include "rstkCoefficientsWriter.h"
#include <itkVectorImageToImageAdaptor.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkDisplacementFieldTransform.h>
#include <itkResampleImageFilter.h>
#include <itkWarpImageFilter.h>

#include "ACWERegistrationMethod.h"
#include "SparseMatrixTransform.h"
#include "DisplacementFieldFileWriter.h"
#include "DisplacementFieldComponentsFileWriter.h"
#include "LevelObserver.h"


using namespace rstk;


namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

const static unsigned int DIMENSION = 3;

typedef float                                                ChannelPixelType;
typedef itk::Image<ChannelPixelType, DIMENSION>              ChannelType;
typedef itk::VectorImage<ChannelPixelType, DIMENSION>        ImageType;
typedef typename ImageType::PixelType                        VectorPixelType;

typedef float                                                ScalarType;
typedef BSplineSparseMatrixTransform<ScalarType,
		                                     DIMENSION, 3u>  TransformType;
typedef TransformType::Pointer                               TransformPointer;
typedef typename TransformType::CoefficientsImageType        CoefficientsType;
typedef typename TransformType::AltCoeffType                 AltCoeffType;

typedef SparseMatrixTransform<ScalarType, DIMENSION >        BaseTransformType;
typedef BaseTransformType::Pointer                           BaseTransformPointer;
typedef typename BaseTransformType::AltCoeffType             AltCoeffType;


typedef itk::ThresholdImageFilter< ChannelType >             ThresholdFilter;
typedef typename ThresholdFilter::Pointer                    ThresholdPointer;

typedef itk::ResampleImageFilter
		         < ChannelType, ChannelType >                ResampleFilter;
typedef typename ResampleFilter::Pointer                     ResamplePointer;
typedef itk::BSplineInterpolateImageFunction
		                      < ChannelType >                DefaultInterpolator;

typedef ACWERegistrationMethod< ImageType, TransformType, ScalarType >   RegistrationType;
typedef typename RegistrationType::Pointer                   RegistrationPointer;
typedef typename RegistrationType::PriorsList                ContourList;
typedef typename RegistrationType::OptimizerType             OptimizerType;
typedef typename RegistrationType::FunctionalType            FunctionalType;
typedef typename FunctionalType::ProbabilityMapType          ProbabilityMapType;
typedef typename OptimizerType::FieldType                    FieldType;
typedef itk::ImageFileReader<ChannelType>                    ImageReader;
typedef itk::ImageFileWriter<ChannelType>                    ImageWriter;
typedef rstk::DisplacementFieldComponentsFileWriter
		                                         <FieldType> ComponentsWriter;
typedef itk::ImageFileWriter< ProbabilityMapType >           ProbabilityMapWriter;
typedef itk::ImageFileWriter< CoefficientsType >             CoefficientsWriter;
typedef rstk::DisplacementFieldFileWriter< FieldType >       FieldWriter;
typedef rstk::CoefficientsWriter< AltCoeffType >             CoeffWriter;

typedef itk::WarpImageFilter
		         < ChannelType, ChannelType, FieldType >     WarpFilter;
typedef typename WarpFilter::Pointer                         WarpFilterPointer;

typedef LevelObserver< RegistrationType >                    LevelObserverType;
typedef typename LevelObserverType::Pointer                  LevelObserverPointer;

typedef typename RegistrationType::PriorsType                PriorType;
typedef rstk::VTKPolyDataWriter<PriorType>                   WriterType;
// typedef itk::MeshFileWriter<PriorType>                       WriterType;

#ifndef NDEBUG
	const static size_t DEFAULT_VERBOSITY = 5;
#else
	const static size_t DEFAULT_VERBOSITY = 1;
#endif

int main(int argc, char *argv[]);

#endif /* REGSEG_H_ */
