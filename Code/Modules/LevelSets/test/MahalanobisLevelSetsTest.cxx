// --------------------------------------------------------------------------
// File:             MahalanobisLevelSetsTest.cxx
// Date:             30/10/2012
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
#include "DisplacementFieldFileWriter.h"

using namespace rstk;

int main(int argc, char *argv[]) {
	typedef itk::Vector<float, 1u>               VectorPixelType;
	typedef itk::Image<VectorPixelType, 3u>      ImageType;
	typedef MahalanobisLevelSets<ImageType>      LevelSetsType;

	typedef LevelSetsType::ContourType     ContourType;
	typedef ContourType::Pointer           ContourDisplacementFieldPointer;
	typedef LevelSetsType::PointType                 PointType;
	typedef LevelSetsType::MeanType                   MeanType;
	typedef LevelSetsType::CovarianceType             CovarianceType;
	typedef LevelSetsType::DeformationFieldType       DeformationFieldType;

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


	LevelSetsType::Pointer ls = LevelSetsType::New();
	ls->SetReferenceImage( im );

	//typename LevelSetsType::ParametersType params;
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


