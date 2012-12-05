// --------------------------------------------------------------------------
// File:             Model1Registration.cxx
// Date:             12/04/2012
// Author:           code@oscaresteban.es (Oscar Esteban, OE)
// Version:          0.1
// License:          BSD
// --------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
// 
// This file is part of ACWEReg
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// * Neither the names of the LTS5-EFPL and the BIT-UPM, nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#ifndef DATA_DIR
#define DATA_DIR "./Data/ModelGeneration/Model1/"
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
#include "MahalanobisLevelSets.h"
#include "GradientDescentLevelSetsOptimizer.h"

using namespace rstk;

int main(int argc, char *argv[]) {
	typedef itk::Vector<float, 1u>                               VectorPixelType;
	typedef itk::Image<VectorPixelType, 3u>                      ImageType;
	typedef MahalanobisLevelSets<ImageType>                      LevelSetsType;

	typedef LevelSetsType::ContourDeformationType                ContourDeformationType;
	typedef ContourDeformationType::Pointer                      ContourDisplacementFieldPointer;
	typedef LevelSetsType::VectorType                            VectorType;
	typedef LevelSetsType::MeanType                              MeanType;
	typedef LevelSetsType::CovarianceType                        CovarianceType;
	typedef LevelSetsType::DeformationFieldType                  DeformationFieldType;

	typedef GradientDescentLevelSetsOptimizer< LevelSetsType >   Optimizer;
	typedef typename Optimizer::Pointer                          OptimizerPointer;

	typedef itk::VTKPolyDataReader< ContourDeformationType >     ReaderType;
	typedef itk::VTKPolyDataWriter< ContourDeformationType >     WriterType;
	typedef itk::ImageFileReader<ImageType>                      ImageReader;
	typedef itk::ImageFileWriter<ImageType>                      ImageWriter;

	typedef itk::VectorImageToImageAdaptor<double,3u>            VectorToImage;

	ImageReader::Pointer r = ImageReader::New();
	r->SetFileName( std::string( DATA_DIR ) + "dtifit__FA.nii.gz" );
	r->Update();
	ImageType::Pointer im = r->GetOutput();

	ReaderType::Pointer polyDataReader = ReaderType::New();
	polyDataReader->SetFileName( std::string( DATA_DIR ) + "ellipse3D-small.vtk" );
	polyDataReader->Update();
	ContourDisplacementFieldPointer initialContour = polyDataReader->GetOutput();

	ReaderType::Pointer polyDataReader2 = ReaderType::New();
	polyDataReader2->SetFileName( std::string( DATA_DIR ) + "wm_prior.vtk" );
	polyDataReader2->Update();
	ContourDisplacementFieldPointer refContour = polyDataReader2->GetOutput();

	typename ContourDeformationType::PointsContainerPointer points = initialContour->GetPoints();
	typename ContourDeformationType::PointsContainerPointer points2 = refContour->GetPoints();
	typename ContourDeformationType::PointsContainerIterator u_it = points->Begin();
	typename ContourDeformationType::PointsContainerIterator u_it2 = points2->Begin();

	VectorType zero = itk::NumericTraits<VectorType>::Zero;
	double initialDistance = 0;
	while( u_it != points->End() ) {
		initialContour->SetPointData( u_it.Index(),zero);
		VectorType dist = u_it.Value() - u_it2.Value();
		initialDistance+= dist.GetSquaredNorm();
		++u_it;
		++u_it2;
	}

	std::cout << "Initial MSE="<< initialDistance / initialContour->GetNumberOfPoints() << std::endl;
	assert( initialContour->GetNumberOfPoints() == refContour->GetNumberOfPoints() );


	MeanType mean1; mean1.Fill(127);
	MeanType mean2; mean2.Fill(255);
	CovarianceType cov; cov.SetIdentity();

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


	LevelSetsType::Pointer ls = LevelSetsType::New();
	ls->SetReferenceImage( im );
	ls->SetContourDeformation( initialContour );
	ls->SetParameters(mean2,cov, true);
	ls->SetParameters(mean1,cov, false);

	OptimizerPointer opt = Optimizer::New();
	opt->SetLevelSetsFunction( ls );
	opt->SetDeformationField( df );
	opt->SetNumberOfIterations(10);
	opt->Start();

	typename ContourDeformationType::ConstPointer resultContour = ls->GetContourDeformation();
	typename ContourDeformationType::PointsContainerConstPointer points3 = resultContour->GetPoints();
	typename ContourDeformationType::PointsContainerConstIterator u_it3 = points3->Begin();
	u_it2 = points2->Begin();

	double finalDistance = 0;
	while( u_it3 != points3->End() ) {
		VectorType dist = u_it3.Value() - u_it2.Value();
		finalDistance+= dist.GetSquaredNorm();
		++u_it2;
		++u_it3;
	}

	std::cout << "Final MSE=" << finalDistance / refContour->GetNumberOfPoints() << std::endl;

	WriterType::Pointer polyDataWriter = WriterType::New();
	polyDataWriter->SetInput( ls->GetContourDeformation() );
	polyDataWriter->SetFileName( "result-registered.vtk" );
	polyDataWriter->Update();




}
