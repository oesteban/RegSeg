/*
 * warp_image.h
 *
 *  Created on: Mar 27, 2014
 *      Author: oesteban
 */

#ifndef WARP_IMAGE_H_
#define WARP_IMAGE_H_


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
#include <itkWarpImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>

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

typedef itk::Vector< ScalarType, DIMENSION >                 VectorType;
typedef itk::Image< VectorType, DIMENSION >                  DisplacementFieldType;
typedef typename DisplacementFieldType::Pointer              DisplacementFieldPointer;
typedef itk::ImageFileReader<DisplacementFieldType>          DisplacementFieldReaderType;
typedef typename DisplacementFieldReaderType::Pointer        DisplacementFieldReaderPointer;

typedef itk::BSplineInterpolateImageFunction< ChannelType, ScalarType >    BSplineInterpolateImageFunction;

typedef itk::WarpImageFilter
		         < ChannelType, ChannelType, DisplacementFieldType >     WarpFilter;
typedef typename WarpFilter::Pointer                         WarpFilterPointer;

int main(int argc, char *argv[]);

#endif /* WARP_IMAGE_H_ */
