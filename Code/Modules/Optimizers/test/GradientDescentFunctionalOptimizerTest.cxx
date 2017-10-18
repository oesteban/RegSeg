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
#include "GradientDescentFunctionalOptimizer.h"

using namespace rstk;

int main(int argc, char *argv[]) {
	typedef itk::Vector<float, 1u>                               VectorPixelType;
	typedef itk::Image<VectorPixelType, 3u>                      ImageType;
	typedef MahalanobisFunctional<ImageType>                      FunctionalType;

	typedef FunctionalType::ContourDeformationType                ContourDeformationType;
	typedef ContourDeformationType::Pointer                      ContourDisplacementFieldPointer;
	typedef FunctionalType::PointType                             PointType;
	typedef FunctionalType::MeanType                              MeanType;
	typedef FunctionalType::CovarianceType                        CovarianceType;
	typedef FunctionalType::DeformationFieldType                  DeformationFieldType;

	typedef GradientDescentFunctionalOptimizer< FunctionalType >   Optimizer;
	typedef typename Optimizer::Pointer                          OptimizerPointer;

	typedef itk::VTKPolyDataReader< ContourDeformationType >     ReaderType;
	typedef itk::VTKPolyDataWriter< ContourDeformationType >     WriterType;
	typedef itk::ImageFileReader<ImageType>                      ImageReader;
	typedef itk::ImageFileWriter<ImageType>                      ImageWriter;

	typedef itk::VectorImageToImageAdaptor<double,3u>            VectorToImage;

	ImageReader::Pointer r = ImageReader::New();
	r->SetFileName( std::string( TEST_DATA_DIR ) + "ellipse3D.nii.gz" );
	r->Update();
	ImageType::Pointer im = r->GetOutput();

	ReaderType::Pointer polyDataReader = ReaderType::New();
	polyDataReader->SetFileName( std::string( TEST_DATA_DIR ) + "ellipse3D-small.vtk" );
	polyDataReader->Update();
	ContourDisplacementFieldPointer initialContour = polyDataReader->GetOutput();

	ReaderType::Pointer polyDataReader2 = ReaderType::New();
	polyDataReader2->SetFileName( std::string( TEST_DATA_DIR ) + "ellipse3D.vtk" );
	polyDataReader2->Update();
	ContourDisplacementFieldPointer refContour = polyDataReader2->GetOutput();

	typename ContourDeformationType::PointsContainerPointer points = initialContour->GetPoints();
	typename ContourDeformationType::PointsContainerPointer points2 = refContour->GetPoints();
	typename ContourDeformationType::PointsContainerIterator u_it = points->Begin();
	typename ContourDeformationType::PointsContainerIterator u_it2 = points2->Begin();

	PointType zero = itk::NumericTraits<PointType>::Zero;
	double initialDistance = 0;
	while( u_it != points->End() ) {
		initialContour->SetPointData( u_it.Index(),zero);
		PointType dist = u_it.Value() - u_it2.Value();
		initialDistance+= dist.GetSquaredNorm();
		++u_it;
		++u_it2;
	}

	std::cout << "Initial MSE="<< initialDistance / initialContour->GetNumberOfPoints() << std::endl;
	assert( initialContour->GetNumberOfPoints() == refContour->GetNumberOfPoints() );

//	CovarianceType cov; cov.SetIdentity();
//	typename FunctionalType::ParametersType params;
//	params.mean[0] = 255;
//	params.mean[1] = 127;
//	params.iCovariance[0] = cov;
//	params.iCovariance[2] = cov;

	DeformationFieldType::Pointer df = DeformationFieldType::New();
	DeformationFieldType::SizeType imSize = im->GetLargestPossibleRegion().GetSize();
	DeformationFieldType::SizeType size; size.Fill(16);
	DeformationFieldType::SpacingType spacing = im->GetSpacing();
	for( size_t i = 0; i<3; i++) spacing[i] = 1.0*imSize[i]/size[i];
	df->SetRegions( size );
	df->SetSpacing( spacing );
	df->SetDirection( im->GetDirection() );
	df->Allocate();
	df->FillBuffer( itk::NumericTraits<DeformationFieldType::PixelType>::Zero );
	std::cout << "Number Of Parameters=" << df->GetLargestPossibleRegion().GetNumberOfPixels() << std::endl;


	FunctionalType::Pointer ls = FunctionalType::New();
	ls->SetReferenceImage( im );
	ls->AddShapePrior( initialContour );

	OptimizerPointer opt = Optimizer::New();
	opt->SetFunctionalFunction( ls );
	//opt->SetDeformationField( df );
	opt->SetNumberOfIterations(10);
	opt->Start();


	typename ContourDeformationType::Pointer resultContour = ls->GetCurrentContourPosition()[0];
	typename ContourDeformationType::PointsContainerConstPointer points3 = resultContour->GetPoints();
	typename ContourDeformationType::PointsContainerConstIterator u_it3 = points3->Begin();
	u_it2 = points2->Begin();

	double finalDistance = 0;
	while( u_it3 != points3->End() ) {
		PointType dist = u_it3.Value() - u_it2.Value();
		finalDistance+= dist.GetSquaredNorm();
		++u_it2;
		++u_it3;
	}

	std::cout << "Final MSE=" << finalDistance / refContour->GetNumberOfPoints() << std::endl;

	WriterType::Pointer polyDataWriter = WriterType::New();
	polyDataWriter->SetInput( ls->GetCurrentContourPosition()[0] );
	polyDataWriter->SetFileName( "result-registered.vtk" );
	polyDataWriter->Update();




}
