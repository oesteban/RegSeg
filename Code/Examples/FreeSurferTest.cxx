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
