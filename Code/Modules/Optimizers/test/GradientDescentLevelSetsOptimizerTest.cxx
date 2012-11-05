// --------------------------------------------------------------------------
// File:             GradientDescentLevelSetsOptimizerTest.cxx
// Date:             01/11/2012
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
	r->SetFileName( std::string( TEST_DATA_DIR ) + "ellipse.nii.gz" );
	r->Update();
	ImageType::Pointer im = r->GetOutput();

	ReaderType::Pointer polyDataReader = ReaderType::New();
	polyDataReader->SetFileName( std::string( TEST_DATA_DIR ) + "ellipse.vtk" );
	polyDataReader->Update();
	ContourDisplacementFieldPointer ellipse = polyDataReader->GetOutput();

	typename ContourDeformationType::PointsContainerPointer points = ellipse->GetPoints();
	typename ContourDeformationType::PointsContainerIterator u_it = points->Begin();

	VectorType zero = itk::NumericTraits<VectorType>::Zero;
	while( u_it != points->End() ) {
		ellipse->SetPointData( u_it.Index(),zero);
		++u_it;
	}

	MeanType mean1; mean1.Fill(127);
	MeanType mean2; mean2.Fill(255);
	CovarianceType cov; cov.SetIdentity();

	DeformationFieldType::Pointer df = DeformationFieldType::New();
	DeformationFieldType::SizeType size = im->GetLargestPossibleRegion().GetSize();
	double factor = 1/3.2;
	for( size_t i=0;i<3;i++) size[i]=(unsigned int) (size[i] * factor);
	DeformationFieldType::SpacingType spacing = im->GetSpacing();
	spacing/=factor;
	df->SetRegions( size );
	df->SetSpacing( spacing );
	df->SetDirection( im->GetDirection() );
	df->Allocate();
	df->FillBuffer( itk::NumericTraits<DeformationFieldType::PixelType>::Zero );
	std::cout << "Number Of Parameters=" << df->GetLargestPossibleRegion().GetNumberOfPixels() << std::endl;


	LevelSetsType::Pointer ls = LevelSetsType::New();
	ls->SetReferenceImage( im );
	ls->SetContourDeformation( ellipse );
	ls->SetParameters(mean2,cov, true);
	ls->SetParameters(mean1,cov, false);

	OptimizerPointer opt = Optimizer::New();
	opt->SetLevelSets( ls );
	opt->Start();



	typedef itk::Image<float,4u> FieldType;
	FieldType::Pointer out = FieldType::New();
	FieldType::SizeType outSize;
	outSize[0] = df->GetLargestPossibleRegion().GetSize()[0];
	outSize[1] = df->GetLargestPossibleRegion().GetSize()[1];
	outSize[2] = df->GetLargestPossibleRegion().GetSize()[2];
	outSize[3] = 3;
	out->SetRegions( outSize );
	FieldType::SpacingType outSpacing;
	outSpacing[0] = df->GetSpacing()[0];
	outSpacing[1] = df->GetSpacing()[1];
	outSpacing[2] = df->GetSpacing()[2];
	outSpacing[3] = 1.0;
	out->SetSpacing( outSpacing );
	out->Allocate();
	out->FillBuffer(0.0);

	float* buffer = out->GetBufferPointer();
	DeformationFieldType::PixelType* vectBuffer = df->GetBufferPointer();
	size_t nPix = df->GetLargestPossibleRegion().GetNumberOfPixels();

	for(size_t pix = 0; pix< nPix; pix++) {
		DeformationFieldType::PixelType val = (*vectBuffer)[pix];
		for( size_t i=0; i<3; i++) {
			unsigned int idx = (i*nPix)+pix;
			buffer[idx] = (float) val[i];
		}
	}
	itk::ImageFileWriter<FieldType>::Pointer w = itk::ImageFileWriter<FieldType>::New();
	w->SetInput( out );
	w->SetFileName( std::string( TEST_DATA_DIR ) + "speed.nii.gz" );
	w->Update();
}
