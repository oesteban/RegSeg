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
