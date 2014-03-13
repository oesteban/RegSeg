// --------------------------------------------------------------------------------------
// File:          bspline_field.h
// Date:          Mar 5, 2014
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

#ifndef BSPLINE_FIELD_H_
#define BSPLINE_FIELD_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkRandomImageSource.h>
#include <itkGaussianImageSource.h>
#include <itkImageRandomNonRepeatingIteratorWithIndex.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkMultiplyImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include "BSplineSparseMatrixTransform.h"
#include "DisplacementFieldComponentsFileWriter.h"

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

static const unsigned int DIMENSION = 3;
static const unsigned int NUM_THREADS = 4;

typedef itk::Image<float, DIMENSION>                         ChannelType;
typedef typename ChannelType::Pointer                        ChannelPointer;
typedef itk::ImageFileReader<ChannelType>                    ReaderType;
typedef typename ReaderType::Pointer                         ReaderPointer;
typedef itk::ImageFileWriter<ChannelType>                    WriterType;
typedef typename WriterType::Pointer                         WriterPointer;


typedef double                                               ScalarType;
typedef rstk::BSplineSparseMatrixTransform
		                          < ScalarType, DIMENSION>   Transform;
typedef Transform::Pointer                                   TPointer;
typedef typename Transform::CoefficientsImageType            CoefficientsType;
typedef typename Transform::CoefficientsImageArray           CoefficientsImageArray;
typedef typename Transform::FieldType                        FieldType;

typedef itk::ImageRandomNonRepeatingIteratorWithIndex
		                             <CoefficientsType>      RandomIterator;
typedef itk::ImageFileWriter< CoefficientsType >             CoefficientsWriterType;
typedef CoefficientsWriterType::Pointer                      CoefficientsWriterPointer;


typedef itk::SmoothingRecursiveGaussianImageFilter< CoefficientsType >
	 	 	 	 	 	 	 	 	 	 	 	 	 	 SmoothingFilterType;
typedef typename SmoothingFilterType::Pointer			 SmoothingFilterPointer;
typedef typename SmoothingFilterType::SigmaArrayType     SigmaArrayType;

typedef itk::MinimumMaximumImageCalculator< CoefficientsType > MaxCalc;
typedef typename MaxCalc::Pointer                              MaxCalcPointer;

typedef itk::StatisticsImageFilter< CoefficientsType >         MedianCalc;
typedef typename MedianCalc::Pointer                           MedianCalcPointer;

typedef itk::MultiplyImageFilter< CoefficientsType >           MultiplyFilter;
typedef typename MultiplyFilter::Pointer                       MultiplyPointer;

typedef itk::SubtractImageFilter< CoefficientsType >           SubtractFilter;
typedef typename SubtractFilter::Pointer                       SubtractPointer;

typedef itk::ResampleImageFilter< ChannelType, ChannelType, ScalarType >    ResampleFilter;
typedef typename ResampleFilter::Pointer                       ResamplePointer;

typedef itk::BSplineInterpolateImageFunction< ChannelType, ScalarType >    BSplineInterpolateImageFunction;

typedef rstk::DisplacementFieldComponentsFileWriter<FieldType> ComponentsWriter;

int main(int argc, char *argv[]);

#endif /* BSPLINE_FIELD_H_ */
