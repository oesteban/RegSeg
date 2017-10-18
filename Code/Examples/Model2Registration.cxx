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
#define DATA_DIR "../Data/ModelGeneration/Model2/"
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
#include <itkOrientImageFilter.h>
#include "MahalanobisFunctional.h"
#include "GradientDescentFunctionalOptimizer.h"

using namespace rstk;

int main(int argc, char *argv[]) {
	typedef itk::Image<itk::Vector<float,1u>, 3u>                ChannelType;
	typedef itk::OrientImageFilter< ChannelType, ChannelType >   OrientFilter;
	typedef MahalanobisFunctional<ChannelType>                    FunctionalType;
	typedef FunctionalType::ContourType                ContourType;
	typedef ContourType::Pointer                      ContourDisplacementFieldPointer;
	typedef FunctionalType::VectorType                            VectorType;
	typedef FunctionalType::MeanType                              MeanType;
	typedef FunctionalType::CovarianceType                        CovarianceType;
	typedef FunctionalType::DeformationFieldType                  DeformationFieldType;

	typedef GradientDescentFunctionalOptimizer< FunctionalType >   Optimizer;
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



	ImageReader::Pointer r = ImageReader::New();
	r->SetFileName( std::string( DATA_DIR ) + "deformed2.nii" );
	r->Update();
	ChannelType::Pointer im = r->GetOutput();
	ChannelType::DirectionType dir = im->GetDirection();

	//ChannelType::DirectionType ident; ident.SetIdentity();
	//OrientFilter::Pointer orient = OrientFilter::New();
	//orient->UseImageDirectionOn();
	//orient->SetDesiredCoordinateDirection( ident );
	//orient->SetInput( im );
	//orient->Update();
	//im = orient->GetOutput();
	//ChannelType::PointType newOrig; newOrig.Fill(0.0);
	//im->SetOrigin( newOrig );

	ImageWriter::Pointer ww = ImageWriter::New();
	ww->SetInput( im );
	ww->SetFileName( "refImage.nii" );
	ww->Update();


	VectorType zero = itk::NumericTraits<VectorType>::Zero;

	ReaderType::Pointer polyDataReader1 = ReaderType::New();
	//polyDataReader1->SetFileName( std::string( DATA_DIR ) + "gmwm_surface.vtk" );
	polyDataReader1->SetFileName( std::string( DATA_DIR ) + "new_surf2_wm.vtk" );
	polyDataReader1->Update();
	ContourDisplacementFieldPointer initialContour1 = polyDataReader1->GetOutput();

	typename ContourType::PointsContainerPointer points1 = initialContour1->GetPoints();
	typename ContourType::PointsContainerIterator u_it1 = points1->Begin();

	while( u_it1 != points1->End() ) {
		initialContour1->SetPointData( u_it1.Index(),zero);
		++u_it1;
	}

	ReaderType::Pointer polyDataReader2 = ReaderType::New();
	//polyDataReader2->SetFileName( std::string( DATA_DIR ) + "pial_surface.vtk" );
	polyDataReader2->SetFileName( std::string( DATA_DIR ) + "new_surf2_gm.vtk" );
	polyDataReader2->Update();
	ContourDisplacementFieldPointer initialContour2 = polyDataReader2->GetOutput();

	typename ContourType::PointsContainerPointer points2 = initialContour2->GetPoints();
	typename ContourType::PointsContainerIterator u_it2 = points2->Begin();

	while( u_it2 != points2->End() ) {
		initialContour2->SetPointData( u_it2.Index(),zero);
		++u_it2;
	}

	// Initialize tissue signatures
	MeanType mean0; // This is WM
	mean0[0] = 0.8;
	CovarianceType cov0;
	cov0(0,0) =  0.003;
	MeanType mean1; // This is GM
	mean1[0] = 0.5;
	CovarianceType cov1;
	cov1(0,0) =  0.004;
	MeanType mean3; // This is CSF
	mean3[0] = 0.21;
	CovarianceType cov3;
	cov3(0,0) =  0.0022;

	typename FunctionalType::ParametersType params1;
	params1.mean[0] = mean0;
	params1.mean[1] = mean1;
	params1.iCovariance[0] = cov0;
	params1.iCovariance[1] = cov1;

	typename FunctionalType::ParametersType params2;
	params2.mean[0] = mean1;
	params2.mean[1] = mean3;
	params2.iCovariance[0] = cov1;
	params2.iCovariance[1] = cov3;

	// Initialize deformation field
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

	// Initialize LevelSet function
	FunctionalType::Pointer ls = FunctionalType::New();
	ls->SetReferenceImage( im );
	ls->AddShapePrior( initialContour1, params1 );
	ls->AddShapePrior( initialContour2, params2 );

	// Connect Optimizer
	OptimizerPointer opt = Optimizer::New();
	opt->SetFunctionalFunction( ls );
	//opt->SetDeformationField( df );
	//opt->SetNumberOfIterations(500);
	//opt->SetAlpha( 50 );
	//opt->SetBeta( 100 );
	//opt->SetStepSize( 0.001 );
	opt->SetNumberOfIterations(5000);
	//opt->SetA( 100 );
	//opt->SetB( 1 );
	opt->SetStepSize( 0.0005 );


	// Start
	opt->Start();

	// Write final result out
	WriterType::Pointer polyDataWriter = WriterType::New();
	polyDataWriter->SetInput( ls->GetCurrentContourPosition()[0] );
	polyDataWriter->SetFileName( "Model2-wm.vtk" );
	polyDataWriter->Update();

	// Write final result out
	WriterType::Pointer polyDataWriter2 = WriterType::New();
	polyDataWriter2->SetInput( ls->GetCurrentContourPosition()[1] );
	polyDataWriter2->SetFileName( "Model2-pial.vtk" );
	polyDataWriter2->Update();

//	DeformationWriter::Pointer w = DeformationWriter::New();
//	w->SetInput( opt->GetDeformationField() );
//	w->SetFileName( "Model2_field.nii.gz" );
//	w->Update();
//
//	DisplacementResamplerType::Pointer p = DisplacementResamplerType::New();
//	p->SetInput( opt->GetDeformationField() );
//	p->SetOutputOrigin( im->GetOrigin() );
//	p->SetOutputSpacing( im->GetSpacing() );
//	p->SetOutputDirection( im->GetDirection() );
//	p->SetSize( im->GetLargestPossibleRegion().GetSize() );
//	//p->SetInterpolator( InterpolatorFunction::New() );
//	p->Update();
//	DeformationWriter::Pointer w2 = DeformationWriter::New();
//	w2->SetInput( p->GetOutput() );
//	w2->SetFileName( "Model2_fieldHD.nii.gz" );
//	w2->Update();
/*
	DeformationFieldType::Pointer dfield = DeformationFieldType::New();
	TransformType::Pointer tf = TransformType::New();
	tf->SetDisplacementField( dfield );

	ResamplerType::Pointer res = ResamplerType::New();
	res->SetInput( im );
	res->SetReferenceImage( im );
	res->SetTransform( tf );
	res->UseReferenceImageOn();
	res->Update();

	ImageWriter::Pointer w3 = ImageWriter::New();
	w3->SetInput( res->GetOutput() );
	w3->SetFileName( "FAunwarped.nii.gz" );
	w3->Update();*/
}


