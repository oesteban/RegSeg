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

#include <itkMaskImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkThresholdImageFilter.h>

#include <itkDisplacementFieldTransform.h>
#include <itkMesh.h>
#include <itkVTKPolyDataReader.h>
#include "rstkVTKPolyDataWriter.h"

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


typedef itk::Image< unsigned char, DIMENSION >               MaskType;
typedef typename MaskType::Pointer                           MaskPointer;
typedef itk::ImageFileWriter<MaskType>                       MaskWriter;
typedef itk::BinaryThresholdImageFilter
		                           < ChannelType, MaskType > Binarize;

typedef double                                               ScalarType;

typedef itk::Vector< ScalarType, DIMENSION >                 VectorType;
typedef itk::Image< VectorType, DIMENSION >                  DisplacementFieldType;
typedef typename DisplacementFieldType::Pointer              DisplacementFieldPointer;
typedef itk::ImageFileReader<DisplacementFieldType>          DisplacementFieldReaderType;
typedef typename DisplacementFieldReaderType::Pointer        DisplacementFieldReaderPointer;

typedef itk::BSplineInterpolateImageFunction< ChannelType, ScalarType >    BSplineInterpolateImageFunction;
typedef itk::NearestNeighborInterpolateImageFunction< ChannelType, ScalarType >    NearestNeighborInterpolateImageFunction;
typedef itk::MaskImageFilter< ChannelType, MaskType, ChannelType > MaskFilter;


typedef itk::ThresholdImageFilter< ChannelType >               ThresholdFilter;
typedef typename ThresholdFilter::Pointer                      ThresholdPointer;

typedef itk::WarpImageFilter
	    < ChannelType, ChannelType, DisplacementFieldType >    WarpFilter;
typedef typename WarpFilter::Pointer                           WarpFilterPointer;

typedef typename itk::DisplacementFieldTransform
		                           < ScalarType, DIMENSION>    TransformType;
typedef typename TransformType::Pointer                        TransformPointer;

typedef itk::Mesh< float, DIMENSION >                          MeshType;
typedef typename MeshType::Pointer                             MeshPointer;
typedef typename MeshType::PointType                           MeshPointType;
typedef typename MeshType::PointsContainer::Iterator           PointsIterator;
typedef itk::VTKPolyDataReader<MeshType>                       MeshReaderType;
typedef typename MeshReaderType::Pointer                       MeshReaderPointer;
typedef rstk::VTKPolyDataWriter<MeshType>                      MeshWriterType;
typedef typename MeshWriterType::Pointer                       MeshWriterPointer;

int main(int argc, char *argv[]);

#endif /* WARP_IMAGE_H_ */
