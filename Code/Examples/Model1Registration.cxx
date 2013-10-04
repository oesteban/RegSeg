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
#define DATA_DIR "../Data/ModelGeneration/Model1/"
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
#include <itkComposeImageFilter.h>
#include <itkVectorResampleImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkDisplacementFieldTransform.h>
#include <itkResampleImageFilter.h>
#include "MahalanobisFunctional.h"
#include "SpectralGradientDescentOptimizer.h"
#include "SpectralADMMOptimizer.h"
//#include "GradientDescentFunctionalOptimizer.h"
//#include "ALOptimizer.h"

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
	//typedef SpectralADMMOptimizer< FunctionalType >              Optimizer;

	typedef typename Optimizer::Pointer                          OptimizerPointer;

	typedef itk::VTKPolyDataReader< ContourType >     ReaderType;
	typedef itk::VTKPolyDataWriter< ContourType >     WriterType;
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


	// Initialize LevelSet function
	FunctionalType::Pointer ls = FunctionalType::New();


	InputToVectorFilterType::Pointer comb = InputToVectorFilterType::New();

	ImageReader::Pointer r = ImageReader::New();
	r->SetFileName( std::string( DATA_DIR ) + "deformed_FA.nii.gz" );
	r->Update();
	comb->SetInput(0,r->GetOutput());

	ImageReader::Pointer r2 = ImageReader::New();
	r2->SetFileName( std::string( DATA_DIR ) + "deformed_MD.nii.gz" );
	r2->Update();

	ChannelType::DirectionType dir; dir.SetIdentity();
	comb->SetInput(1,r2->GetOutput());
	comb->Update();
	ImageType::Pointer im = comb->GetOutput();
	//im->SetDirection( dir );

	ls->SetReferenceImage( im );

	ReaderType::Pointer polyDataReader0 = ReaderType::New();
	//polyDataReader0->SetFileName( std::string( DATA_DIR ) + "fixed.csf.vtk" );
	polyDataReader0->SetFileName( std::string( DATA_DIR ) + "surf/roi_csf_smoothed.nii_1_converted.vtk" );
	polyDataReader0->Update();
	ls->AddShapePrior( polyDataReader0->GetOutput() );

	//ReaderType::Pointer polyDataReader1 = ReaderType::New();
	//polyDataReader1->SetFileName( std::string( DATA_DIR ) + "fixed.wm.vtk" );
	//polyDataReader1->Update();
	//ls->AddShapePrior( polyDataReader1->GetOutput() );
    //
	//ReaderType::Pointer polyDataReader2 = ReaderType::New();
	//polyDataReader2->SetFileName( std::string( DATA_DIR ) + "fixed.gm.vtk" );
	//polyDataReader2->Update();
	//ls->AddShapePrior( polyDataReader2->GetOutput() );

	// Connect Optimizer
	OptimizerPointer opt = Optimizer::New();
	opt->SetGridSize( 6 );
	opt->SetFunctional( ls );
	opt->SetNumberOfIterations(10);
	opt->SetAlpha( 0.01 );
	opt->SetStepSize( 0.1 );

	// Start
	opt->Start();

	// Write final result out
    size_t nCont = ls->GetCurrentContours().size();
    for ( size_t contid = 0; contid < nCont; contid++) {
        std::stringstream ss;
        ss << "final-cont0" << contid << ".vtk";
    	WriterType::Pointer polyDataWriter = WriterType::New();
    	polyDataWriter->SetInput( ls->GetCurrentContours()[contid] );
    	polyDataWriter->SetFileName( "deformed2-wm.vtk" );
    	polyDataWriter->Update();
    }
}
