// --------------------------------------------------------------------------------------
// File:          FreeSurferTest.cxx
// Date:          Mar 22, 2013
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.5.5
// License:       GPLv3 - 29 June 2007
// --------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of RegSeg
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


#include "MahalanobisFunctional.h"
#include "SpectralGradientDescentOptimizer.h"
#include "SpectralADMMOptimizer.h"

using namespace rstk;

int main(int argc, char *argv[]) {
	typedef itk::Image<float, 3u>                                ChannelType;
	typedef itk::Vector<float, 2u>                               VectorPixelType;
	typedef itk::Image<VectorPixelType, 3u>                      ImageType;
	typedef itk::ComposeImageFilter< ChannelType,ImageType >     InputToVectorFilterType;

	typedef MahalanobisFunctional<ImageType>                      FunctionalType;
	typedef FunctionalType::ContourType                           ContourType;
	typedef ContourType::Pointer                                 ContourDisplacementFieldPointer;
	typedef FunctionalType::MeanType                              MeanType;
	typedef FunctionalType::CovarianceType                        CovarianceType;
	typedef FunctionalType::FieldType                             DeformationFieldType;

	typedef SpectralGradientDescentOptimizer< FunctionalType >   Optimizer;
	typedef typename Optimizer::Pointer                          OptimizerPointer;

	typedef itk::MeshFileReader< ContourType >                   ReaderType;
	typedef itk::MeshFileWriter< ContourType >                   WriterType;
	typedef itk::ImageFileReader<ChannelType>                    ImageReader;
	typedef itk::ImageFileWriter<ChannelType>                    ImageWriter;
	typedef itk::ImageFileWriter<DeformationFieldType>           DeformationWriter;

	typedef itk::VectorImageToImageAdaptor<double,3u>            VectorToImage;


	typedef itk::Image< itk::Vector<double,3u>, 3u >             NewField;
	typedef itk::VectorResampleImageFilter
			<DeformationFieldType,NewField, double >             DisplacementResamplerType;
	typedef itk::BSplineInterpolateImageFunction
			                <DeformationFieldType>               InterpolatorFunction;

	typedef itk::DisplacementFieldTransform<double, 3u>          TransformType;

	typedef itk::ResampleImageFilter<ChannelType,ChannelType,double>    ResamplerType;




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

	// Initialize LevelSet function
	FunctionalType::Pointer ls = FunctionalType::New();
	ls->SetReferenceImage( comb->GetOutput() );
	//ls->AddShapePrior( initialContour2 );
	//ls->AddShapePrior( initialContour3 );

	// Connect Optimizer
	OptimizerPointer opt = Optimizer::New();
	opt->SetFunctional( ls );
	opt->SetNumberOfIterations(15);

	std::string fnames[3] = { "init.csf.vtk", "init.lh.white.vtk", "init.lh.pial.vtk" };

	for ( size_t n = 0; n<3; n++ ) {
		ReaderType::Pointer polyDataReader0 = ReaderType::New();
		polyDataReader0->SetFileName( std::string( DATA_DIR ) + fnames[n] );
		polyDataReader0->Update();
		ContourDisplacementFieldPointer initialContour0 = polyDataReader0->GetOutput();
		ls->AddShapePrior( initialContour0 );
	}

	// Start
	opt->Start();

	// Write final result out
	WriterType::Pointer polyDataWriter = WriterType::New();
	polyDataWriter->SetInput( ls->GetCurrentContours()[0] );
	polyDataWriter->SetFileName( "registered.white.vtk" );
	polyDataWriter->Update();

	DisplacementResamplerType::Pointer disResampler = DisplacementResamplerType::New();
	disResampler->SetInput( opt->GetCurrentDisplacementField() );
	disResampler->SetOutputOrigin(    im->GetOrigin()  );
	disResampler->SetOutputSpacing(   im->GetSpacing() );
	disResampler->SetSize(            im->GetLargestPossibleRegion().GetSize() );
	disResampler->SetOutputDirection( im->GetDirection() );
	disResampler->Update();

	//DeformationWriter::Pointer w = DeformationWriter::New();
	//w->SetInput( disResampler->GetOutput() );
	//w->SetFileName( "estimated_field.nii.gz" );
	//w->Update();

	TransformType::Pointer t = TransformType::New();
	t->SetDisplacementField( disResampler->GetOutput() );


	typename ResamplerType::Pointer res = ResamplerType::New();
	res->SetInput( im );
	res->SetTransform( t );
	res->SetOutputParametersFromImage( im );
	res->Update();

	ImageWriter::Pointer ww = ImageWriter::New();
	ww->SetInput( res->GetOutput() );
	ww->SetFileName( "FA_resampled.nii.gz" );
	ww->Update();




}
