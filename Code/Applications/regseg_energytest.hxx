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

#ifndef SOURCE_DIRECTORY__APPLICATIONS_REGSEG_ENERGYTEST_HXX_
#define SOURCE_DIRECTORY__APPLICATIONS_REGSEG_ENERGYTEST_HXX_

#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include <itkArray.h>
#include <itkVector.h>
#include <itkVectorImage.h>
#include <itkQuadEdgeMesh.h>
#include <itkVTKPolyDataReader.h>
#include <itkImageFileReader.h>
#include <itkOrientImageFilter.h>
#include <itkComposeImageFilter.h>

#include "DownsampleAveragingFilter.h"
#include "MultilabelBinarizeMeshFilter.h"
#include "rstkVTKPolyDataWriter.h"
#include "ComponentsFileWriter.h"
#include "EnergyCalculatorFilter.h"
#include "MahalanobisDistanceModel.h"

#include <jsoncpp/json/json.h>

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

const static unsigned int Dimension = 3u;

typedef double                                                    MeasureType;
typedef float                                                     ChannelPixelType;
typedef itk::Image<ChannelPixelType, Dimension>                   ChannelType;
typedef itk::VectorImage<ChannelPixelType, Dimension>             ReferenceImageType;
typedef typename ReferenceImageType::PixelType                    VectorPixelType;
typedef itk::ComposeImageFilter<ChannelType, ReferenceImageType>  InputToVectorFilterType;

typedef float                                                     PointValueType;
typedef itk::Vector< PointValueType, Dimension >                  VectorType;
typedef itk::QuadEdgeMesh< VectorType, Dimension >                VectorContourType;
typedef itk::VectorImage< PointValueType, Dimension >             ProbmapType;

typedef typename ReferenceImageType::Pointer                      ReferencePointer;
typedef typename ReferenceImageType::DirectionType                DirectionType;
typedef typename ReferenceImageType::SizeType                     SizeType;
typedef typename ReferenceImageType::PointType                    PointType;
typedef typename ReferenceImageType::SpacingType                  SpacingType;
typedef itk::ContinuousIndex<PointValueType, Dimension>           ContinuousIndex;


typedef rstk::MultilabelBinarizeMeshFilter< VectorContourType >   BinarizeMeshFilterType;
typedef typename BinarizeMeshFilterType::Pointer                  BinarizeMeshFilterPointer;
typedef typename BinarizeMeshFilterType::OutputImageType          BinarizationImageType;
typedef typename BinarizeMeshFilterType::InputMeshContainer       InputMeshContainer;
typedef typename BinarizeMeshFilterType::OutputComponentType      SegmentationType;

typedef rstk::DownsampleAveragingFilter
		< BinarizationImageType, ProbmapType >                    DownsampleFilter;
typedef typename DownsampleFilter::Pointer                        DownsamplePointer;

typedef itk::VTKPolyDataReader< VectorContourType >               ReaderType;
typedef rstk::VTKPolyDataWriter< VectorContourType >              WriterType;
typedef itk::ImageFileReader<ChannelType>                         ImageReader;
typedef itk::ImageFileWriter<SegmentationType>                    SegmentationWriter;
typedef rstk::ComponentsFileWriter<ProbmapType>                   ImageWriter;

typedef itk::OrientImageFilter< ReferenceImageType, ReferenceImageType >  Orienter;
typedef itk::OrientImageFilter< SegmentationType, SegmentationType >      SegmentationOrienter;
typedef itk::OrientImageFilter< ProbmapType, ProbmapType >                ProbmapsOrienter;
typedef itk::OrientImageFilter< ChannelType, ChannelType >                ChannelOrienter;

typedef rstk::MahalanobisDistanceModel< ReferenceImageType >      EnergyModelType;
typedef typename EnergyModelType::Pointer                         EnergyModelPointer;

typedef rstk::EnergyCalculatorFilter<ReferenceImageType>          EnergyFilter;
typedef typename EnergyFilter::Pointer                            EnergyFilterPointer;
typedef typename EnergyFilter::PriorsImageType                    PriorsImageType;
typedef itk::Array< MeasureType >                                 MeasureArray;
typedef typename PriorsImageType::Pointer                         PriorsImagePointer;
typedef typename PriorsImageType::PixelType                       PriorsPixelType;
typedef typename PriorsImageType::InternalPixelType               PriorsValueType;




int main(int argc, char *argv[]);


#endif /* SOURCE_DIRECTORY__APPLICATIONS_REGSEG_ENERGYTEST_HXX_ */
