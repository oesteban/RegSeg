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

#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR "./"
// data source http://code.google.com/p/v3dsolver/source/browse/data/
#endif

#include <itkVector.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkQuadEdgeMesh.h>
#include <itkVTKPolyDataReader.h>
#include <itkVTKPolyDataWriter.h>
#include <itkVectorImageToImageAdaptor.h>
#include "MahalanobisFunctional.h"
#include "DisplacementFieldFileWriter.h"

using namespace rstk;

int main(int argc, char *argv[]) {
	typedef itk::Vector<float, 1u>               VectorPixelType;
	typedef itk::Image<VectorPixelType, 3u>      ImageType;
	typedef MahalanobisFunctional<ImageType>      FunctionalType;

	typedef FunctionalType::ContourType     ContourType;
	typedef ContourType::Pointer           ContourDisplacementFieldPointer;
	typedef FunctionalType::PointType                 PointType;
	typedef FunctionalType::MeanType                   MeanType;
	typedef FunctionalType::CovarianceType             CovarianceType;
	typedef FunctionalType::DeformationFieldType       DeformationFieldType;

	typedef itk::VTKPolyDataReader< ContourType >     ReaderType;
	typedef itk::VTKPolyDataWriter< ContourType >     WriterType;
	typedef itk::ImageFileReader<ImageType>                      ImageReader;
	typedef itk::ImageFileWriter<ImageType>                      ImageWriter;
	typedef rstk::DisplacementFieldFileWriter<DeformationFieldType> Writer;

	typedef itk::VectorImageToImageAdaptor<double,3u>            VectorToImage;

	ImageReader::Pointer r = ImageReader::New();
	r->SetFileName( std::string( TEST_DATA_DIR ) + "ellipse.nii.gz" );
	r->Update();
	ImageType::Pointer im = r->GetOutput();

	ReaderType::Pointer polyDataReader = ReaderType::New();
	polyDataReader->SetFileName( std::string( TEST_DATA_DIR ) + "ellipse.vtk" );
	polyDataReader->Update();
	ContourDisplacementFieldPointer ellipse = polyDataReader->GetOutput();

	typename ContourType::PointsContainerPointer points = ellipse->GetPoints();
	typename ContourType::PointsContainerIterator u_it = points->Begin();

	PointType zero = itk::NumericTraits<PointType>::Zero;
	while( u_it != points->End() ) {
		ellipse->SetPointData( u_it.Index(),zero);
		++u_it;
	}

	MeanType mean1; mean1.Fill(127);
	MeanType mean2; mean2.Fill(255);
	CovarianceType cov; cov.SetIdentity();

	DeformationFieldType::Pointer df = DeformationFieldType::New();
	DeformationFieldType::SizeType size = im->GetLargestPossibleRegion().GetSize();
	double factor = 0.25;
	for( size_t i=0;i<3;i++) size[i]=(unsigned int) (size[i] * factor);
	DeformationFieldType::SpacingType spacing = im->GetSpacing();
	spacing=spacing/factor;
	df->SetRegions( size );
	df->SetSpacing( spacing );
	df->SetDirection( im->GetDirection() );
	df->Allocate();
	df->FillBuffer( itk::NumericTraits<DeformationFieldType::PixelType>::Zero );
	std::cout << "Number Of Parameters=" << df->GetLargestPossibleRegion().GetNumberOfPixels() << std::endl;


	FunctionalType::Pointer ls = FunctionalType::New();
	ls->SetReferenceImage( im );

	//typename FunctionalType::ParametersType params;
	//params.mean[0] = 255;
	//params.mean[1] = 127;
	//params.cov = cov;
	//params.iCovariance[2] = cov;
	ls->AddShapePrior( ellipse );

	Writer::Pointer writer = Writer::New();
	writer->SetFileName( std::string( TEST_DATA_DIR ) + "shapegradientmap.nii.gz" );
	writer->SetInput( ls->GetShapeGradients() );
	writer->Update();
}


