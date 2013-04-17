// --------------------------------------------------------------------------------------
// File:          FreeSurferTest.cxx
// Date:          Mar 22, 2013
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
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


/*
 * Test 1 - Read vtk mesh, save it as FreeSurfer binary mesh and open it with tksurfer
 * Test 2 - Read previous FreeSurfer mesh, launch registration and save result as FreeSurfer binary mesh.
 */
#ifndef DATA_DIR
#define DATA_DIR "../Data/FreeSurferTest/"
#endif

#include <itkVector.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkQuadEdgeMesh.h>
#include <itkMeshFileReader.h>
#include <itkMeshFileWriter.h>
#include <itkVectorImageToImageAdaptor.h>
#include <itkComposeImageFilter.h>
#include <itkVectorResampleImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkDisplacementFieldTransform.h>
#include <itkResampleImageFilter.h>
#include "MahalanobisLevelSets.h"
#include "GradientDescentLevelSetsOptimizer.h"

using namespace rstk;

int main(int argc, char *argv[]) {
	typedef itk::Image<float, 3u>                                ChannelType;
	typedef itk::Vector<float, 2u>                               VectorPixelType;
	typedef itk::Image<VectorPixelType, 3u>                      ImageType;
	typedef itk::ComposeImageFilter< ChannelType,ImageType >     InputToVectorFilterType;

	typedef MahalanobisLevelSets<ImageType>                      LevelSetsType;
	typedef LevelSetsType::ContourDeformationType                ContourDeformationType;
	typedef ContourDeformationType::Pointer                      ContourDisplacementFieldPointer;
	typedef LevelSetsType::VectorType                            VectorType;
	typedef LevelSetsType::MeanType                              MeanType;
	typedef LevelSetsType::CovarianceType                        CovarianceType;
	typedef LevelSetsType::DeformationFieldType                  DeformationFieldType;

	typedef GradientDescentLevelSetsOptimizer< LevelSetsType >   Optimizer;
	typedef typename Optimizer::Pointer                          OptimizerPointer;

	typedef itk::MeshFileReader< ContourDeformationType >     ReaderType;
	typedef itk::MeshFileWriter< ContourDeformationType >     WriterType;
	typedef itk::ImageFileReader<ChannelType>                      ImageReader;
	typedef itk::ImageFileWriter<ChannelType>                      ImageWriter;
	typedef itk::ImageFileWriter<DeformationFieldType>           DeformationWriter;

	typedef itk::VectorImageToImageAdaptor<double,3u>            VectorToImage;
	typedef itk::VectorResampleImageFilter
			<DeformationFieldType,DeformationFieldType,double>   DisplacementResamplerType;
	typedef itk::BSplineInterpolateImageFunction
			                <DeformationFieldType>               InterpolatorFunction;
	typedef itk::DisplacementFieldTransform<float, 3u>           TransformType;

	typedef itk::ResampleImageFilter<ChannelType,ChannelType,float>    ResamplerType;




	InputToVectorFilterType::Pointer comb = InputToVectorFilterType::New();

	ImageReader::Pointer r = ImageReader::New();
	r->SetFileName( std::string( DATA_DIR ) + "FA.nii.gz" );
	r->Update();
	comb->SetInput(0,r->GetOutput());

	ImageReader::Pointer r2 = ImageReader::New();
	r2->SetFileName( std::string( DATA_DIR ) + "MD.nii.gz" );
	r2->Update();
	comb->SetInput(1,r2->GetOutput());
	comb->Update();

	ChannelType::Pointer im = r->GetOutput();


	ReaderType::Pointer polyDataReader = ReaderType::New();
	polyDataReader->SetFileName( std::string( DATA_DIR ) + "lh.white.vtk" );
	polyDataReader->Update();
	ContourDisplacementFieldPointer initialContour = polyDataReader->GetOutput();

	typename ContourDeformationType::PointsContainerPointer points = initialContour->GetPoints();
	typename ContourDeformationType::PointsContainerIterator u_it = points->Begin();

	VectorType zero = itk::NumericTraits<VectorType>::Zero;
	while( u_it != points->End() ) {
		initialContour->SetPointData( u_it.Index(),zero);
		++u_it;
	}

	// Initialize tissue signatures
	MeanType mean1; // This is GM
	mean1[0] = 0.11941234;
	mean1[1] = 0.00089523;
	CovarianceType cov1;
	cov1(0,0) =  5.90117156e-04;
	cov1(0,1) = -1.43226633e-06;
	cov1(1,0) = -1.43226633e-06;
	cov1(1,1) =  1.03718252e-08;
	MeanType mean2; // This is WM
	mean2[0] = 7.77335644e-01;
	mean2[1] = 6.94673450e-04;
	CovarianceType cov2;
	cov2(0,0) =  4.85065832e-03;
	cov2(0,1) = -6.89610616e-06;
	cov2(1,0) = -6.89610616e-06;
	cov2(1,1) =  1.02706528e-08;

	typename LevelSetsType::ParametersType params;
	params.mean[0] = mean2;
	params.mean[1] = mean1;
	params.iCovariance[0] = cov2;
	params.iCovariance[2] = cov1;

	// Initialize LevelSet function
	LevelSetsType::Pointer ls = LevelSetsType::New();
	ls->SetReferenceImage( comb->GetOutput() );
	ls->AddShapePrior( initialContour, params );

	// Connect Optimizer
	OptimizerPointer opt = Optimizer::New();
	opt->SetLevelSetsFunction( ls );
	opt->SetNumberOfIterations(15);

	// Start
	opt->Start();

	// Write final result out
	WriterType::Pointer polyDataWriter = WriterType::New();
	polyDataWriter->SetInput( ls->GetCurrentContourPosition()[0] );
	polyDataWriter->SetFileName( "registered.white.vtk" );
	polyDataWriter->Update();

	DeformationWriter::Pointer w = DeformationWriter::New();
	w->SetInput( opt->GetDeformationField() );
	w->SetFileName( "estimated_field.nii.gz" );
	w->Update();
}
