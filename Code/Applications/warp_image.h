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
#include <itkResampleImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>

#include <itkMaskImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkThresholdImageFilter.h>

#include <itkDisplacementFieldTransform.h>
#include <itkMesh.h>
#include <itkVTKPolyDataReader.h>
#include "rstkVTKPolyDataWriter.h"
#include "CompositeMatrixTransform.h"
#include "BSplineSparseMatrixTransform.h"

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


typedef float                                                ScalarType;
typedef rstk::BSplineSparseMatrixTransform
		                          < ScalarType, DIMENSION>   BSplineTransform;
typedef BSplineTransform::Pointer                            BSplineTransformPointer;
typedef rstk::CompositeMatrixTransform
		                          < ScalarType, DIMENSION>   CompositeTransform;
typedef CompositeTransform::Pointer                          CompositeTransformPointer;


typedef typename BSplineTransform::CoefficientsImageType     CoefficientsType;
typedef typename BSplineTransform::CoefficientsImageArray    CoefficientsImageArray;
typedef typename BSplineTransform::FieldType                 FieldType;
typedef typename BSplineTransform::VectorType                VectorType;

typedef typename FieldType::Pointer              			 DisplacementFieldPointer;
typedef typename FieldType::ConstPointer              		 DisplacementFieldConstPointer;
typedef itk::ImageFileReader<FieldType>          			 DisplacementFieldReaderType;
typedef typename DisplacementFieldReaderType::Pointer        DisplacementFieldReaderPointer;

typedef itk::MaskImageFilter< ChannelType, MaskType, ChannelType > MaskFilter;

typedef itk::ResampleImageFilter< ChannelType, ChannelType, ScalarType >    ResampleFilter;
typedef typename ResampleFilter::Pointer                       ResamplePointer;
typedef itk::BSplineInterpolateImageFunction< ChannelType, ScalarType >    ResampleInterpolator;

typedef itk::ResampleImageFilter
		                 < MaskType, MaskType, ScalarType >    ResampleMaskFilter;
typedef typename ResampleMaskFilter::Pointer                   ResampleMaskPointer;
typedef itk::NearestNeighborInterpolateImageFunction< MaskType, ScalarType >    ResampleMaskInterpolator;

typedef itk::WarpImageFilter
	    < MaskType, MaskType, FieldType >                      WarpMaskFilter;
typedef typename WarpMaskFilter::Pointer                       WarpMaskFilterPointer;
typedef itk::NearestNeighborInterpolateImageFunction< MaskType >    WarpMaskInterpolator;

typedef itk::WarpImageFilter
	    < ChannelType, ChannelType, FieldType >                WarpFilter;
typedef typename WarpFilter::Pointer                           WarpFilterPointer;
typedef itk::BSplineInterpolateImageFunction< ChannelType >    WarpInterpolator;

typedef itk::ThresholdImageFilter< ChannelType >               ThresholdFilter;
typedef typename ThresholdFilter::Pointer                      ThresholdPointer;

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
